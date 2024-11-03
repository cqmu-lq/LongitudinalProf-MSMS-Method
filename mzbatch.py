import os,subprocess,pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyopenms import MSExperiment,MzMLFile,PeakPickerHiRes
import pymzml

# 本文件为公共方法，分为以下模块：文件读写，数据检查，数据生成，数据计算，数据提取

"""=============================
文件读写
"""
# 文件复制
import shutil
def copyfile(infile,outfile):
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        makedir(outdir)
        while True:
            if os.path.exists(outdir):
                break
    shutil.copy(infile,outfile)
    
# 创建文件夹：文件夹不存在则创建
def makedir(path):
    if not os.path.exists(path):
        subprocess.Popen('mkdir "{}"'.format(path), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
def make_folder(folder_path):
    if not os.path.exists(folder_path):
            os.makedirs(folder_path)

# 文件路径是否存在: 不存在创建，创建好了才返回True
def dir_is_exist(filepath):
    if os.path.exists(filepath):
        return True
    else:
        os.makedirs(filepath)
        while True:
            if os.path.exists(filepath):
                break
        return True

# 写入pickle文件
def to_pkl(data,w_path):
    with open(w_path,'wb') as f:
        pickle.dump(data,f)
# 读取pickle文件
def read_pkl(r_path):
    with open(r_path, 'rb') as f:
        data=pickle.load(f)
    return data
        
# 读取mzML数据
def readmzML(filepath):
    exp=MSExperiment()
    MzMLFile().load(filepath, exp)
    return exp

# 解析mzML文件
def get_TIC_from_File(file_path,tic_pkl):
    '''
    return [ {'rt','scanID',,,,},    ]
    '''
    if os.path.exists(tic_pkl):
        TIC = read_pkl(tic_pkl)
    else:
        # print(type(file_path))
        # 载入文件
        if isinstance(file_path,str):
            exp=MSExperiment()
            MzMLFile().load(file_path, exp)
        elif isinstance(file_path,MSExperiment):
            exp=file_path
        else:
            return

        profiled_exp=MSExperiment()
        temp_cen_list=[]
        for spec in exp:
            if spec.getType()==2:
                profiled_exp.addSpectrum(spec)
            else:
                temp_cen_list.append(spec)
        # 归一化数据
        centroided_exp = MSExperiment()
        if profiled_exp.size()>0:
            PeakPickerHiRes().pickExperiment(exp, centroided_exp)

        # print("centroided_exp.size()  ",centroided_exp.size())
        # print( "len(temp_cen_list) ",len(temp_cen_list))
        if len(temp_cen_list)>0:
            for spec in temp_cen_list:
                centroided_exp.addSpectrum(spec)
        centroided_exp.sortSpectra()
        # print("centroided_exp.size()  ",centroided_exp.size())
        
        # pymzml读，获取isolation window target m/z
        TIC2 = pymzml.run.Reader(file_path)

        TIC=[]
        for idx,spec in enumerate( centroided_exp):
            spectrum={}
            spectrum['ID']=idx+1
            spectrum['RT']=spec.getRT()
            spectrum['Polarity']=spec.getInstrumentSettings().getPolarity()
            spectrum['MSlevel']=spec.getMSLevel()
            cur_id = idx if spectrum['MSlevel']==1 else cur_id
            cur_RT = spectrum['RT'] if spectrum['MSlevel']==1 else cur_RT
            spectrum['MS1_id']=cur_id if spectrum['MSlevel']==2 else None
            spectrum['MS1_RT']=cur_RT if spectrum['MSlevel']==2 else None
            # spectrum['pecursor']= None if spectrum['MSlevel']==1 else spec.getPrecursors()[0].getMZ()
            
            # 就修改这一句
            spectrum['pecursor']= None if spectrum['MSlevel']==1 else TIC2[idx+1]['MS:1000827'] 
            
            # spectrum['up_offset']= None if spectrum['MSlevel']==1 else spec.getPrecursors()[0].getIsolationWindowUpperOffset()
            # spectrum['low_offset']= None if spectrum['MSlevel']==1 else spec.getPrecursors()[0].getIsolationWindowLowerOffset()
            spectrum['peaks']=list(zip(*spec.get_peaks()))
            spectrum['TIC']=np.sum(np.array(spectrum['peaks'])[:,1])
            TIC.append(spectrum)
        # print(list(zip(*centroided_exp[1].get_peaks())))
        to_pkl(TIC,tic_pkl)
    return TIC

# 合并文件夹下的所有csv结果
def marge_data_fun(folder_path,qianzhui,begin_0,end_0,step_0,outfile):
    all_files_list = ["{}/{}_{}_{}.csv".format(folder_path,qianzhui,i,i+step_0) for i in range(begin_0,end_0,step_0)]
    merged_data_df = pd.DataFrame()
    for file in all_files_list:
        df = pd.read_csv(file,index_col=0)
        merged_data_df = pd.concat([merged_data_df, df], ignore_index=True)
    print(merged_data_df.shape, outfile)
    merged_data_df.to_excel(outfile,encoding='utf8')
"""=============================
数据检查
"""

# 数据校验：1. 统计一二级谱图数量；2. 检查是不是每个一级谱图后面都有4个二级谱图
def check_data(tic_list):
    # 统计一二级谱图数量
    MS1_count,MS2_count=0,0
    for tic_ in tic_list:
        if tic_['MSlevel']==1:
            MS1_count+=1
        else:
            MS2_count+=1
    count_all=len(tic_list)
    print('1. 总数、一级和二级数量分别为：',count_all,MS1_count,MS2_count)
    i=0
    id_list = []
    for one_tic in tic_list:
        if one_tic['MSlevel']==1:
            id=one_tic['ID']
            end=id+4 if id+4<count_all else count_all

            if not [tic_list[id_]['MSlevel'] for id_ in range(id,end)]==[2,2,2,2]:
                id_list.append(id-1)
    print("2. 输出id后面不是跟4个的index（已经减1）：",id_list)

"""=============================
数据生成
"""
# 将26386按照 ‘文件名-能量-pos/neg.mzML’的格式生成文件列表--（仅使用于第一次），第二次文件名为 “{}-NCE{}-{}”
def get_filename_list(filepre_list):
    filename_list = []
    ev_list = ['1', '10', '100', '15', '20', '25', '30', '35', '40', '50', '60','70', '80']
    for filepre in filepre_list:
        for ev in ev_list:
            for tail in ['pos','neg']:
                filename_list.append('{}-{}-{}'.format(filepre,ev,tail))
    return filename_list

# 说明：生成一级图文件
# 输入：temp_df 小印给的化合物列表，tic_MS1_df是tic_list中一级谱图列表
# pos_neg是正负离子，MS1_2dict_file是一级谱图数据的二维字典，key是cas，内容是rt和intensity
def make_MS1_data(tic_list,temp_df,tic_MS1_df,pos_neg,MS1_2dict_file):
    if os.path.exists(MS1_2dict_file):
        MS1_2dict = read_pkl(MS1_2dict_file)
    else:
        MS1_2dict={}
        rt_list = []
        # 先构建好字典         
        for temp_index in temp_df.index:
            cas = temp_df.loc[temp_index,'CAS No.']
            MS1_2dict[cas] = {'rt':rt_list,'intensity':[]}
        # 往字典中放数据
        for tic_MS1_index in tic_MS1_df.index:
            tic_one_dict = tic_list[tic_MS1_index]
            d_df= pd.DataFrame(tic_one_dict['peaks'],columns=['mass','intensity'])
            rt_list.append(tic_one_dict['RT'])
            for i in temp_df.index:
                cas = temp_df.loc[i,'CAS No.']
                # 旧
                # mz = temp_df.loc[i,'m/z-'+pos_neg]
                # 240409修改
                mz = temp_df.loc[i,'m/z']
                
                filter_df = d_df[d_df['mass'].apply(lambda x: less_5ppm_true(abs(cal_ppm(mz,x))))]
                if filter_df.shape[0] != 0:
                    MS1_2dict[cas]['intensity'].append(filter_df.sum()[1])
                else:
                    MS1_2dict[cas]['intensity'].append(0)
        to_pkl(MS1_2dict,MS1_2dict_file)
    return MS1_2dict

"""=============================
数据计算
5ppm改10ppm，只需要改这里就行，改5个数
"""
# 计算ppm
def cal_ppm(expected_value,real_value):
    return (expected_value-real_value)/expected_value*1000000

# 是否在 ±5ppm 以内
# 说明：一级计算还用的是5ppm，二级用的10ppm，为了最大限度减少代码修改，增加less_5ppm_true方法，用于一级5ppm
def less_5ppm_true(ppm_value):
    return True if abs(ppm_value)<=5 else False
def less_5ppm(ppm_value):
    return True if abs(ppm_value)<=10 else False
def less_10ppm(ppm_value):
    return True if abs(ppm_value)<=10 else False

# 计算理论mz上下限
def get_mz_uplimit(expected_value,diff=10):
    return expected_value+diff*expected_value/1000000
def get_mz_lowlimit(expected_value,diff=10):
    return expected_value-diff*expected_value/1000000
def get_mz_range(expected_value,diff=10):
    mz_up = get_mz_uplimit(expected_value,diff)
    mz_low = get_mz_lowlimit(expected_value,diff)
    return mz_low,mz_up
# 矩阵加速
def get_mz_range_matrix(expected_values, diff=10):
    # Convert expected_values to a numpy array
    expected_values = np.array(expected_values)
    
    # Calculate the mz range using matrix operations
    mz_low = expected_values - diff * expected_values / 1000000
    mz_up = expected_values + diff * expected_values / 1000000
    
    # Combine the lower and upper limits into a list of tuples
    mz_range_list = np.vstack((mz_low, mz_up)).T
    return mz_range_list.tolist()

# 判断mz是否大于上下限，返回bool型的array数组，用于取代apply实现加速
def mz_is_less5(data_Series,theo_mz):
    mz_low,mz_up = get_mz_range(theo_mz)
    return (data_Series.values >= mz_low)&(data_Series.values <= mz_up)
# 判断df的某一列，在上下限之间
def between_low_up(data_Series,low_,up_):
    return (data_Series.values >= low_)&(data_Series.values <= up_)
# 判断大于是否大于上限
def mz_bigger_uplimit(expected_value,real_value):
    mz_uplimit = get_mz_uplimit(expected_value)
    return real_value <= mz_uplimit

"""=============================
数据提取
"""

#说明：获取上下限
#先生成阈值，找最小值，如果最小值<1%,取1%，如果最小值>1%,最小值/最大值+1%
def get_threshold(data_df):
    min_ratio = (data_df['intensities_list'].min()-100) /data_df['intensities_list'].max()
    if min_ratio<0.01:
        return 0.01
    else:
        return min_ratio+0.01
# 获取最大intentisy所在的区间
def get_max_rt_limt(limit_2list,max_intensity_rt):
    for low_,up_ in limit_2list:
        if (low_<=max_intensity_rt<=up_):
            return [low_,up_]
            break
            
def get_rt_limit(rt_list,intensities_list,refer_rt=-1,name=""):
    # 一级谱图
    mz_RT_df = pd.DataFrame({'rt_list':rt_list,'intensities_list':intensities_list})

    # 如果一级谱全为0，说明这个raw文件中，没有找到这个化合物
    no0_df = mz_RT_df[mz_RT_df['intensities_list']!=0]
    if no0_df.shape[0]==0:
        explain = 'intensity全是0'
        return [-1,-1],explain
    
    no0_df = no0_df.reset_index(drop=True)
    max_intensity = no0_df['intensities_list'].max()
    max_intensity_rt = no0_df[no0_df['intensities_list'] == max_intensity].iloc[0,0]
    # print(max_intensity_rt)
    
    # 获取阈值
    threshold = get_threshold(no0_df)
    threshold_intensity = max_intensity * threshold
    
    # 得到分段的df
    threshold_df = mz_RT_df['intensities_list']>=threshold_intensity
    # threshold_df = mz_RT_df['intensities_list']>0
    threshold_diff = threshold_df.diff()
    threshold_diff.loc[0]=True # 接入头尾
    threshold_diff.loc[threshold_diff.shape[0]-1]=True
    mz_threshold_diff_df = mz_RT_df[threshold_diff]
    
    # 数据分段，并获得上下限, 合并[(1,2),(3,4)] 和 [(2,3)]
    limit_2list = list(zip(mz_threshold_diff_df.iloc[::2, 0], mz_threshold_diff_df.iloc[1::2, 0]))
    limit_2list.extend(list(zip(mz_threshold_diff_df.iloc[1::2, 0], mz_threshold_diff_df.iloc[2::2, 0])))
    limit_2list = sorted(limit_2list, key=lambda x: x[0])
    # print(limit_2list)
    # print('max_intensity_rt',max_intensity_rt)
    lower_limit_rt,upper_limit_rt = [-1,-1]
    explain=''
    max_rt_low,max_rt_up = get_max_rt_limt(limit_2list,max_intensity_rt)
    for low_,up_ in limit_2list:
        
        # 区间内的df
        range_df = mz_RT_df[(mz_RT_df['rt_list']>=low_) & (mz_RT_df['rt_list']<=up_)]
        # 判断区间是否全部小于最小阈值,不要当前区间。下面还有种情况，卡5e5
        if (range_df.iloc[:-1,1]<threshold_intensity).all():
            # print("全部小于阈值")
            continue
        
        # 有参考rt        
        if (refer_rt!=-1):
            # 峰宽是否超过60:
            if (up_ - low_)>60:
                # refer_rt 在最大intensity±20以内，取最大±20；否则取 refer_rt±20
                if (max_intensity_rt-20<=refer_rt<=max_intensity_rt+20):
                    lower_limit_rt,upper_limit_rt = max_intensity_rt-20,max_intensity_rt+20
                    explain='峰宽 > 60, refer_rt在max±20以内'
                    break
                else:
                    lower_limit_rt,upper_limit_rt = refer_rt-20,refer_rt+20
                    explain='峰宽 > 60, refer_rt在max以外，只选refer_rt的±20'
                    break
            else:
                # 峰宽小于60，找refer_rt和 max_intensity_rt 都在里面的区间；然后找refer_rt所在的区间；最后找最大intensity±12有rt的
                if (low_<=refer_rt<=up_) & (low_<=max_intensity_rt<=up_):
                    lower_limit_rt,upper_limit_rt = low_-6,up_+6
                    explain='峰宽<60, refer_rt和max都在区间，区间±6'
                    break
                elif (low_<=refer_rt<=up_) :
                    # 最大在rt±6以内，用最大intensity所在的区间±6
                    if (low_-6<=max_intensity_rt<=up_+6):
                        lower_limit_rt,upper_limit_rt = max_rt_low-6,max_rt_up+6
                        explain='峰宽<60, 只有refer_rt在区间,max在±6，取max区间的±6'
                        break
                    # 当前区间的最大强度大于5e5
                    range_max_intensity = range_df['intensities_list'].max()
                    if range_max_intensity>=5e5:
                        # range_max_rt = range_df[range_df['intensities_list'] == range_max_intensity].iloc[0,0]
                        lower_limit_rt,upper_limit_rt = low_-6,up_+6
                        explain='峰宽<60, 只有refer_rt在区间,max不在±6,range_max>5e5，取当前区间±6'
                        break
                elif (low_<=max_intensity_rt<=up_) :
                    if (low_-6<refer_rt<up_+6):
                        lower_limit_rt,upper_limit_rt = low_-6,up_+6
                        explain='峰宽<60, max过滤，refer_rt在其±6以内，区间±6'
                        break
                    
                
        # 没有rt，直接从6-30中找
        # if (refer_rt==-1) & (low_-12<=max_intensity_rt<=up_+12):
        if (refer_rt==-1):
            lower_limit_rt,upper_limit_rt = 6,30
            explain='没有refer_rt，直接找6-30中间的'
            break
                
    if ([lower_limit_rt,upper_limit_rt] == [-1,-1]):
        # 3.23增加规则
        if max_intensity>=1e6:
            lower_limit_rt,upper_limit_rt = max_rt_low-6,max_rt_up+6
            explain='没有找到区间，max大于1e6，直接max所在区间±6'
        else:
            explain='彻底没找到'
            return [-1,-1],explain

    return [lower_limit_rt,upper_limit_rt],explain

#=================
#===最优谱图判断===
#================

# 说明：获取当前质谱数据中的母离子峰，如果没有返回[-1,-1]
# 输入：质谱数据，分子基本信息
# 输出：[母离子质量,强度]
def get_real_val(intensity_df,theo_mz):
    # temp_df = intensity_df[intensity_df['Mass'].apply(lambda x:less_5ppm(cal_ppm(theo_mz,x)))]
    temp_df = intensity_df[mz_is_less5(intensity_df['Mass'],theo_mz)]
    
    if temp_df.empty:
        return [-1,-1]
    else:
        temp_max=temp_df['Intensity'].max()
        temp_max_df = temp_df[temp_df['Intensity']==temp_max]
        return [temp_max_df.iloc[0,0],temp_max_df.iloc[0,1]]
# 说明：获取基峰。返回[-1,-1]说明剔除上限峰以后数据空了
# 输入：质谱数据，分子基本信息
# 输出：[基峰质量,强度]
def get_basepeak(data_df,theo_mz=""):
    try:
        if theo_mz != "":
            # 计算母离子峰的上限
            mz_uplimit = get_mz_uplimit(theo_mz)
            # 剔除超过上限的峰
            data_df = data_df[data_df['Mass'].values<=mz_uplimit]
        temp_max=data_df['Intensity'].max()
        temp_df=data_df[data_df['Intensity']==temp_max]
        return [temp_df.iloc[0,0],temp_df.iloc[0,1]]
    except Exception as e:
        print('get_basepeak error:',e)
        return [-1,-1]
# 说明：判断基峰 Mass 是否合理
# 基峰和母离子峰差值在4-13；19-25的不合理。 3.11日 改为(4-13,20-25)。12.13日 改 [3,13],[21,24]
# 输出：[是否合理，差值]
def basepeak_is_reasonable(molpeak_mass,basepeak_mass):
    mass_dif=molpeak_mass-basepeak_mass
    if(mass_dif<-1)|(mass_dif<=0.99999)|(3<=mass_dif<=13)|(21<=mass_dif<=24):
        return [False,mass_dif]
    else:
        return [True,mass_dif]
# 说明：判断母离子峰是否大于基峰
# 输入：母峰强度，基峰强度，阈值5%
# 输出：[是否大于5%，百分比]
def greater_base_5per(molpeak_intensity,basepeak_intnesity, Threshold):
    percent=round(molpeak_intensity/basepeak_intnesity,4)*100
    if percent >=Threshold:
        return True,percent
    else:
        return False,percent
# 说明：计算碎片峰数量，强度大于指定基峰5%的峰
# 输入：读取的每个质谱文件（两列，mass和Intensity），碎片峰的判断阈值
# 输出：大于指定阈值的峰合集df
def get_fragmentpeak_df(intensity_df,basepeak_intensity,Threshold,theo_mz):
    fragmentpeak_df=intensity_df[(intensity_df['Intensity']/basepeak_intensity)*100>=Threshold]
    # print('===')
    # print(fragmentpeak_df)
    # 去除不合理峰，用basepeak_is_reasonable方法
    fragmentpeak_df = fragmentpeak_df[fragmentpeak_df['Mass'].apply(lambda x: basepeak_is_reasonable(theo_mz,x)).apply(lambda x: pd.Series(x))[0]]
    return fragmentpeak_df

# 说明：碎片峰是否大于等于1
def fragmentpeak_num_greater1(fragmentpeak_df):
    return True if fragmentpeak_df.shape[0]>1 else False

# 说明：获取输入文件夹下 MS2的文件名
def get_ms2_file(foldername):
    # 先判断文件夹是否存在
    if not os.path.exists(foldername):
        return -1
    ms2_filename_list = []
    for i in os.listdir(foldername):
        if '_MS2_' in i:
            ms2_filename_list.append(i)
    if ms2_filename_list == []:
        return -1
    else:
        return ms2_filename_list
    
# 获取当前二级图中基峰的intensity，注意要删除大于母离子上限的峰
# 输入文件路径
def get_ms2_base(filepath,theo_mz):
    if not os.path.exists(filepath):
        return 0
    try:
        ms2_data = pd.read_csv(filepath,index_col=0)
        # 计算母离子峰的上限
        mz_uplimit = get_mz_uplimit(theo_mz)
        # 剔除超过上限的峰
        ms2_data = ms2_data[ms2_data['Mass']<=mz_uplimit]
        # return mzbatch.get_real_val(ms2_data,theo_mz)[1]
    except Exception:
        print(filepath)
        return 0
    return get_basepeak(ms2_data)[1]

"""=============================
相似度计算
"""
from sklearn.metrics.pairwise import cosine_similarity
def round_df(input_df,rate=10):
    if input_df.shape[0]>0:
        input_df['Mass']=round(input_df['Mass']*rate)
        max_indexes = input_df.groupby(by=['Mass'])['Intensity'].idxmax()
        input_df = input_df.loc[max_indexes]
        input_df['Mass']=input_df['Mass'].astype(int)
        # input_df.set_index(['Mass'],inplace=True)
        return input_df
    else:
        input_df['Mass']=input_df['Mass'].astype(int)
        return input_df

def layer_round_df(input_df):
    input_low_df = input_df[input_df['Mass']<=100]
    input_high_df = input_df[input_df['Mass']>100]

    input_low_df = round_df(input_low_df,1000)
    input_high_df = round_df(input_high_df,100)
    input_high_df = round_df(input_high_df,10)

    input_df = pd.concat([input_low_df,input_high_df])
    # print(input_low_df)

    return input_df
# rate这个参数已经没用了
def layer_read_data(filepath,mz,rate=10):
    # data_df = pd.read_csv(filepath,index_col=0)
    data_df = pd.read_csv(filepath,usecols=["Mass","Intensity"])
    # 计算母离子峰的上限，并剔除超过上限的峰
    mz_uplimit = get_mz_uplimit(mz)
    data_df = data_df[data_df['Mass']<=mz_uplimit]
    
    data_df = layer_round_df(data_df)
    # a_df.columns=['Mass','5300-03-8_35_pos1']
    data_df.set_index(['Mass'],inplace=True)
    return data_df

# rate这个参数已经没用了
def read_data(filepath,mz,rate=100):
    # data_df = pd.read_csv(filepath,index_col=0)
    data_df = pd.read_csv(filepath,usecols=["Mass","Intensity"])
    # 计算母离子峰的上限，并剔除超过上限的峰
    mz_uplimit = get_mz_uplimit(mz)
    data_df = data_df[data_df['Mass']<=mz_uplimit]
    
    data_df = round_df(data_df,rate)
    # a_df.columns=['Mass','5300-03-8_35_pos1']
    data_df.set_index(['Mass'],inplace=True)
    return data_df

def change_data(data_df,mz):
    # 计算母离子峰的上限，并剔除超过上限的峰
    mz_uplimit = get_mz_uplimit(mz)
    data_df = data_df[data_df['Mass']<=mz_uplimit]
    
    data_df = layer_round_df(data_df)
    # a_df.columns=['Mass','5300-03-8_35_pos1']
    data_df.set_index(['Mass'],inplace=True)
    return data_df

#======相似度计算======
def fu_jaccard_similarity_matrix(X, Y=None):
    # If Y is None, use X for both arguments (i.e., compare all pairs within X)
    if Y is None:
        Y = X

    # Ensure input is numpy arrays
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)

    # Compute dot products
    dot_product = np.dot(X, Y.T)
    
    # Compute the sum of squares for each vector
    sum_of_squares_X = np.sum(X ** 2, axis=1)[:, np.newaxis]  # Make column vector
    sum_of_squares_Y = np.sum(Y ** 2, axis=1)  # Make row vector
    
    # Compute the Jaccard similarity matrix
    denominator = sum_of_squares_X + sum_of_squares_Y - dot_product
    jaccard_index = dot_product / denominator
    
    return jaccard_index

import numpy as np
from sklearn.utils import check_array
def entropy_similarity(x, y=None):
    """
    Calculate the entropy similarity between two vectors or matrices.

    Parameters:
    - x: array-like, shape (n_samples, n_features) or (n_features,)
        First input array or vector.
    - y: array-like, shape (n_samples, n_features) or (n_features,), optional
        Second input array or vector. If None, the similarity will be computed
        between the rows of x.

    Returns:
    - Similarity matrix or vector, depending on the shape of the inputs.
    """
    def cal_Ipie(x_one, s_x):
        if s_x >= 3:
            return x_one
        else:
            w = 0.25 + s_x * 0.25
            return np.power(x_one, w)

    def compute_entropy(x):
        x = np.clip(x, 1e-10, None) # 将0都替换为1e-10
        x_sum = np.sum(x)
        x_norm = x / x_sum
        s_x = -np.sum([i * np.log(i) for i in x_norm])
        # print(x_norm, s_x)
        return x_norm, s_x

    def compute_entropy_pie(x, s_x):
        x = np.clip(x, 1e-10, None)
        x_pie = np.array([cal_Ipie(i, s_x) for i in x])
        x_pie_sum = np.sum(x_pie)
        x_pie_norm = x_pie / x_pie_sum
        s_x_pie = -np.sum([i * np.log(i) for i in x_pie_norm])
        return x_pie_norm, s_x_pie

    x = check_array(x, accept_sparse=False, ensure_2d=False)
    if y is not None:
        y = check_array(y, accept_sparse=False, ensure_2d=False)

    if x.ndim == 1:
        # Handle case where x is a vector
        x_norm, s_x = compute_entropy(x)
        if y is None:
            # Calculate similarity with itself
            y_norm, s_y = compute_entropy(x)
            x_pie_norm, s_x_pie = compute_entropy_pie(x, s_x)
            y_pie_norm, s_y_pie = compute_entropy_pie(y, s_y)
            ab_pie_norm = (x_pie_norm + y_pie_norm) / 2
            s_ab_pie = -np.sum([i * np.log(i) for i in ab_pie_norm])
            entropy_sim = 1 - ((2 * s_ab_pie - s_x_pie - s_y_pie) / np.log(4))
            return entropy_sim
        else:
            # Calculate similarity between x and y
            y_norm, s_y = compute_entropy(y)
            x_pie_norm, s_x_pie = compute_entropy_pie(x, s_x)
            y_pie_norm, s_y_pie = compute_entropy_pie(y, s_y)
            ab_pie_norm = (x_pie_norm + y_pie_norm) / 2
            s_ab_pie = -np.sum([i * np.log(i) for i in ab_pie_norm])
            entropy_sim = 1 - ((2 * s_ab_pie - s_x_pie - s_y_pie) / np.log(4))
            return entropy_sim

    elif x.ndim == 2:
        # Handle case where x is a matrix
        n_samples_x = x.shape[0]
        n_samples_y = y.shape[0] if y is not None else n_samples_x

        if y is None:
            similarity_matrix = np.zeros((n_samples_x, n_samples_x))
            for i in range(n_samples_x):
                for j in range(i, n_samples_x):
                    entropy_sim = entropy_similarity(x[i], x[j])
                    similarity_matrix[i, j] = entropy_sim
                    similarity_matrix[j, i] = entropy_sim
            return similarity_matrix
        else:
            similarity_matrix = np.zeros((n_samples_x, n_samples_y))
            for i in range(n_samples_x):
                for j in range(n_samples_y):
                    entropy_sim = entropy_similarity(x[i], y[j])
                    similarity_matrix[i, j] = entropy_sim
            return similarity_matrix

#以下是-1列为未知物#
#240802修改：增加熵相似度计算
def cal_pos_sim(mz_filter_df, sim_method = "cosine"):
    pos_sim_df = mz_filter_df[mz_filter_df[-1]!=0]
    if pos_sim_df.shape[0]==0:
        return [0 for i in range(pos_sim_df.shape[1]-1)]
    if sim_method == "cosine":
        pos_matrix = cosine_similarity(pos_sim_df.T)
    elif sim_method == "jaccard":
        pos_matrix = fu_jaccard_similarity_matrix(pos_sim_df.T)
    elif sim_method == "entropy":
        pos_matrix = entropy_similarity(pos_sim_df.T)
    pos_sim_list = list(np.round(pos_matrix[0][1:],4))  
    return pos_sim_list



#反向模式计算
def cal_neg_sim(standard_ms2_df,mz_filter_df, sim_method = "cosine"):
    standard_columns_list = list(standard_ms2_df.columns)
    neg_sim_list = []
    for temp_index in standard_columns_list:
        temp_df = mz_filter_df.loc[:,[-1,temp_index]]
        temp_df = temp_df[temp_df[temp_index]!=0]
        # 有些二级图是空的,相似度取-1
        if temp_df.shape[0]==0:
            neg_sim_list.append(-1)
        else:
            if sim_method == "cosine":
                neg_sim_list.append(cosine_similarity(temp_df.T)[0][1])
            elif sim_method == "jaccard":
                neg_sim_list.append(fu_jaccard_similarity_matrix(temp_df.T)[0][1])
            elif sim_method == "entropy":
                neg_sim_list.append(entropy_similarity(temp_df.T)[0][1])
    neg_sim_list = list(np.round(neg_sim_list,4))
    
    return neg_sim_list

# 双边计算
def cal_both_sim(mz_filter_df, sim_method = "cosine"):
    both_sim_df = mz_filter_df
    if sim_method == "cosine":
        both_matrix = cosine_similarity(both_sim_df.T)
    elif sim_method == "jaccard":
        both_matrix = fu_jaccard_similarity_matrix(both_sim_df.T)
    elif sim_method == "entropy":
        both_matrix = entropy_similarity(both_sim_df.T)
    both_sim_list = list(np.round(both_matrix[0][1:],4))
    return both_sim_list


# 转为相对丰度
def to_relative_abundance(data_df):
    # 找到每列的最大值
    max_values = data_df.max()
    # 将每一列除以该列的最大值
    result = data_df.div(max_values)
    return result

"""=============================
绘图
"""
# 绘制总离子流图
def draw_tic(exp):
    tic = exp.calculateTIC()
    retention_times, intensities = tic.get_peaks()
    plt.plot(retention_times, intensities)
    plt.title('TIC')
    plt.xlabel('time (s)')
    plt.ylabel('intensity (cps)')
    plt.show()

# 绘图
# 输入：横坐标，纵坐标，标题，图类型（影响横纵轴标签），x标签范围，垂直的红线
def plot_pic(x_list, y_list,title,pic_type,label_lim=None,rt_line_list=None):
    title=str(title)
    plt.plot(x_list, y_list)

    plt.title(title)
    if pic_type=='rt':
        xlabel,ylabel='time (s)','intensity (cps)'
    elif pic_type=='mz':
        xlabel,ylabel='mass','intensity'   
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    if label_lim is not None:
        plt.xlim(*label_lim)
    if rt_line_list is not None:
        for rt_lint in rt_line_list:
            plt.axvline(rt_lint,color='red')
    plt.show()
    
"""=============================
不用方法
"""

# 说明：（旧），根据生成一级谱，这个数据是没有归一化的，该方法不用了获取指定mz 5ppm以内的强度和保留时间
# 输入：mzML文件读取内容exp，指定mz
# 输出：一级谱图的绘制数据
def getRT_intensity(exp,mz,MS1_file):
    if os.path.exists(MS1_file):
        mz_RT_df=pd.read_csv(MS1_file,index_col=0)
        rt_list,intensities_list=mz_RT_df['rt_list'],mz_RT_df['intensities_list']
    else:
        rt_list= []
        intensities_list=[]
        for spec in exp:
            if spec.getMSLevel() == 1:
                mass_array,intensity_array=spec.get_peaks()
                d_df=pd.DataFrame({'mass':mass_array,'intensity':intensity_array})
                d_df=d_df[d_df['mass'].apply(lambda x: less_5ppm(abs(cal_ppm(mz,x))))]
                rt_list.append(spec.getRT())
                intensities_list.append(d_df.sum()[1])
        # print(len(rt_list),len(intensities_list))
                # for mass_, intensity_ in zip(*spec.get_peaks()):
                #     if less_5ppm(abs(cal_ppm(mz,mass_))):
                #         intensitie_list.append(intensity_)
        # 保存到文件
        mz_RT_df=pd.DataFrame({'rt_list':rt_list,'intensities_list':intensities_list})
        # mz_RT_df.to_csv(MS1_file)
    return rt_list, intensities_list