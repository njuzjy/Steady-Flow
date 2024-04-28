import pandas as pd
import numpy as np

##根据水位，计算断面的过流面积、湿周以及水利半径
def Cal_Sec(H,Start_Dis,Elevation):
    '''
    :param
    H: 给定的水位，浮点型
    Start_Dis: 起点距,np.array类型
    Elevation: 高程，np.array类型
    :return:
    Area：过水面积
    Wet_Perimeter：湿周
    Water_Width：过水宽度
    R：水利半径
    '''
    ##定义变量
    Area = 0.0           ##过水面积
    R = 0.0              ##水利半径
    Wet_Perimeter = 0.0  ##湿周
    Water_Width = 0.0    ##过水宽度

    Point_diff = np.zeros(Start_Dis.shape)
    Point_diff = np.where(H>=Start_Dis,1,0)
    #当给定的水位H小于河底时，提示错误
    if H <= np.min(Elevation):
        print("水位小于河底")
    for i in range(len(Start_Dis)-1):
        #计算过流面积
        if H > Elevation[i] and H <= Elevation[i+1]:
            Area = Area + (H - Elevation[i])**2 * (Start_Dis[i+1] - Start_Dis[i]) / (Elevation[i+1] - Elevation[i]) / 2     #计算过流面积，单位m2
            Wet_Perimeter = Wet_Perimeter + ((Start_Dis[i+1] - Start_Dis[i]) ** 2 + (Elevation[i+1] - Elevation[i]) ** 2) ** (1/2) * (H - Elevation[i]) / (Elevation[i+1] - Elevation[i])       #计算湿周，单位m
            Water_Width = Water_Width + abs(H - Elevation[i]) * (Start_Dis[i+1] - Start_Dis[i]) / (Elevation[i+1] - Elevation[i])       #计算过水宽度，单位m

        elif H <= Elevation[i] and H > Elevation[i+1]:
            Area = Area + (H - Elevation[i+1])**2 * (Start_Dis[i] - Start_Dis[i+1]) / (Elevation[i+1] - Elevation[i]) / 2   #计算过流面积，单位m2
            Wet_Perimeter = Wet_Perimeter + ((Start_Dis[i+1] - Start_Dis[i]) ** 2 + (Elevation[i+1] - Elevation[i]) ** 2) ** (1/2) * (H - Elevation[i+1]) / (Elevation[i] - Elevation[i+1])     #计算湿周，单位m
            Water_Width = Water_Width + (H - Elevation[i+1]) * (Start_Dis[i] - Start_Dis[i+1]) / (Elevation[i+1] - Elevation[i])      #计算过水宽度，单位m

        elif H >= Elevation[i] and H >= Elevation[i+1]:
            Area = Area + (2*H - Elevation[i] - Elevation[i+1]) * (Start_Dis[i+1] - Start_Dis[i]) /2                        #计算过流面积，单位m2
            Wet_Perimeter = Wet_Perimeter + ((Start_Dis[i+1] - Start_Dis[i]) ** 2 + (Elevation[i+1] - Elevation[i]) ** 2) ** (1/2)                                                              #计算湿周，单位m
            Water_Width = Water_Width + Start_Dis[i+1] - Start_Dis[i]      #计算过水宽度，单位m

        else:
            Area = Area                                                                                                     #计算过流面积，单位m2
            Wet_Perimeter = Wet_Perimeter                                                                                                                                                       #计算湿周，单位m
            Water_Width = Water_Width

    ##计算水利半径R
    R = Area / Water_Width

    return (Area,Wet_Perimeter,Water_Width,R)

##恒定流计算核心程序一
def Cal_Steady_Flow(Hu,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L):
    '''
    :param
    Hu: 上游断面水位
    Hd: 下游断面水位
    Sec_u: 上游断面数据，包含起点距，高程
    Sec_d: 下游断面数据，包含起点距，高程
    n_u: 上游断面糙率
    n_d: 下游断面糙率
    Alpha_u: 上游断面修正系数
    Beta_u: 上游断面的局部阻水系数
    Alpha_d: 下游断面修正系数
    Beta_d: 下游断面的局部阻水系数
    L：两个断面之间的距离
    :return:
    '''
    g = 9.81
    #获取上下游断面数据
    (Area_u,_,_,R_u) = Cal_Sec(Hu,Sec_u[0],Sec_u[1])
    (Area_d,_,_,R_d) = Cal_Sec(Hd,Sec_d[0],Sec_d[1])

    #计算上下游流速
    V_u = 0.0 ; V_d = 0.0
    V_u = Q_u / Area_u ; V_d = Q_d / Area_d

    #计算局部水头损失Hfik
    Hj = 0.0
    Hj = (Beta_u + Beta_d) / 2 * (V_d ** 2 / (2 * g) - V_u ** 2 / (2 * g))

    # 计算流量模数K
    K_u = 0.0; K_d = 0.0
    K_u =  1 / n_u * R_u **(1 / 6) * Area_u * R_u ** (1 / 2)
    K_d =  1 / n_d * R_d **(1 / 6) * Area_d * R_d ** (1 / 2)

    # 计算沿程水头损失
    Hf = 0.0
    Hf = (Q_u ** 2 + Q_d ** 2) / 2 * L * (1/ (K_u ** 2) + 1 / (K_d ** 2)) / 2

    # 计算主函数
    err = Hu + Alpha_u * V_u ** 2 / (2 * g) - Hd - Alpha_d * V_d ** 2 / (2 * g) - Hf - Hj

    return err

##恒定流计算核心程序二
##本公式的计算来源于论文《城市排水系统的整合与优化研究》
def Steady_Flow(Hu,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L):
    '''
    :param
    Hu: 上游断面水位
    Hd: 下游断面水位
    Sec_u: 上游断面数据，包含起点距，高程
    Sec_d: 下游断面数据，包含起点距，高程
    n_u: 上游断面糙率
    n_d: 下游断面糙率
    Alpha_u: 上游流速系数
    Beta_u: 上游断面的局部阻水系数
    Alpha_d: 下游流速系数
    Beta_d: 下游断面的局部阻水系数
    L：两个断面之间的距离
    :return:
    '''
    pass




## 读取恒定流模型输入文件
def Read_ModelParam(dir):
    '''
    :param
    dir: 文件路径
    :return:
    Params:文件中的参数
    断面名称	流量系数	主槽糙率	局损系数	距离
    '''
    Params = pd.read_table(dir)
    return Params['断面名称'],Params['流量系数'],Params['主槽糙率'],Params['局损系数'],Params['距离']

## 根据断面编号，获取断面的里程及高程信息
def Get_Section_Info(DMMC,dir):
    '''
    DMBH: 断面编号，字符串
    dir: 文件路径
    :return:
    Sec_Info：断面信息，np.array
    '''
    ##定义列表，DMLC用来存储断面里程，DMGC用来存储断面高程
    DMLC = []; DMGC = []
    ##读取txt文件
    with open(dir,encoding="utf8") as f1:
        contents = f1.readlines()
        for index,content in enumerate(contents):
            if content.split('\t')[0] == DMMC:
                DMLC = contents[index+1].split('\t')
                DMGC = contents[index+2].split('\t')
    #print(SKSW,SKKR)
    ##去除掉断面里程和断面高程里面的空值，包括空格，tab键和换行符
    DMLC = [LC for LC in DMLC if LC.strip()]
    DMGC = [GC for GC in DMGC if GC.strip()]
    #print(DMLC,DMGC)
    Sec_Info = Str_Change_Array(DMLC,DMGC)
    return Sec_Info,min(DMGC)

##读取初始条件CSTJ文件
def Read_CSTJ_Params(dir):
    '''
    :param
    dir: 文件路径
    :return:
    SYLL：上游流量
    XYSW：下游水位值
    '''
    Params = pd.read_table(dir)
    SYLL = Params['SYLL']
    XYSW = Params['XYSW']
    #print(SYLL,XYSW)
    return SYLL[0],XYSW[0]

##将Str元素转为float数组
def Str_Change_Array(list1,list2):
    '''
    list1: 字符串型list
    list2: 字符串型list
    return:
    Array：二维的array
    '''
    W_Level = list(map(float,list1))
    Relate_W = list(map(float,list2))
    Array = np.array((W_Level,Relate_W))
    return Array

##二分法求解每个断面的值
def Dichotomy(H_max,H_min,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L):
    '''
    :param H_max:初步迭代计算的最大值
    :param H_min:初步迭代计算的最小值
    :param Hd:
    :param Sec_u:
    :param Sec_d:
    :param n_u:
    :param n_d:
    :param Alpha_u:
    :param Beta_u:
    :param Alpha_d:
    :param Beta_d:
    :param Q_u:
    :param Q_d:
    :param L:
    :return:
    H_Mid：最优结果
    '''
    result_min = Cal_Steady_Flow(H_min,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L)
    result_max = Cal_Steady_Flow(H_max,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L)
    #print(result_min,result_max)
    while H_min <= H_max:
        H_Mid = (H_min + H_max) / 2
        result = Cal_Steady_Flow(H_Mid,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L)
        #print(result)
        if abs(result) <= 0.01:
            break
        if result_min * result < 0:
            H_min = H_min
            H_max = H_Mid
            result_max = Cal_Steady_Flow(H_max,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L)
        if result_max * result < 0:
            H_min = H_Mid
            H_max = H_max
            result_min = Cal_Steady_Flow(H_min,Hd,Sec_u,Sec_d,n_u,n_d,Alpha_u,Beta_u,Alpha_d,Beta_d,Q_u,Q_d,L)
    return H_Mid


##测试程序
if __name__ == '__main__':
    #读取input文件
    Input_dir = "D:\学习\计算模型\恒定非均匀流\分段试算法\模型计算\东源港\input.txt"
    Section_dir = "D:\学习\计算模型\恒定非均匀流\分段试算法\模型计算\东源港\Section.txt"
    CSTJ_dir = "D:\学习\计算模型\恒定非均匀流\分段试算法\模型计算\东源港\CSTJ.txt"
    DMs,LLXSs,Ns,Betas,Distances = Read_ModelParam(Input_dir)
    SYLL,XYSW = Read_CSTJ_Params(CSTJ_dir)
    SW_Output = [XYSW]
    for i in range(len(DMs)-1):
        ##获取下游断面的信息
        Sec_Info_d,Sec_min_d = Get_Section_Info(DMs[i],Section_dir)
        #获取上游断面的信息
        Sec_Info_u,Sec_min_u = Get_Section_Info(DMs[i+1],Section_dir)
        if i == 0:
            #print(float(Sec_min_u)+50.0,float(Sec_min_u)+0.001)
            H_u = Dichotomy(float(Sec_min_u)+50.0,float(Sec_min_u)+0.001,float(XYSW),Sec_Info_u,Sec_Info_d,Ns[i+1],Ns[i],1.5,Betas[i+1],1.5,Betas[i],LLXSs[i+1]*float(SYLL),LLXSs[i]*float(SYLL),Distances[i+1])
            #print(H_u)
            SW_Output.append(H_u)
        else:
            H_u = Dichotomy(float(Sec_min_u)+50.0,float(Sec_min_u)+0.001,H_u,Sec_Info_u,Sec_Info_d,Ns[i+1],Ns[i],1.5,Betas[i+1],1.5,Betas[i],LLXSs[i+1]*float(SYLL),LLXSs[i]*float(SYLL),Distances[i+1])
            #print(H_u)
            SW_Output.append(H_u)


    ##写入out文件
    with open(r'D:\学习\计算模型\恒定非均匀流\分段试算法\模型计算\东源港\out.txt','w') as f:
        f.write("DMBM	SW\n")
        for i in range(len(SW_Output)):
            f.write(str(DMs[i])+str('\t')+str(round(SW_Output[i],3))+str('\n'))


