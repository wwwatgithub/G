# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 16:20:50 2022

@author: EVEN
"""

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from sympy import *
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


N = 400
M = 400
H = 20#发射点与接收点的高度
H1 = 2
wavelength = 0.129 #波长
ds = gdal.Open("E:\\实验数据\\DEM-平原-丘陵-高山\\ASTGTM2_N34E114_dem_Clip_Res1.tif")
#ds = gdal.Open("E:\\arcgis\\arcmap\\dem1_region1_Resample1_Clip1.tif")
rows = ds.RasterXSize
cols = ds.RasterYSize
bands = ds.RasterCount 
band1 = ds.GetRasterBand(1)
band1_array = band1.ReadAsArray() 
ds_array = ds.ReadAsArray()
#裁剪栅格为指定大小
ds_array = ds_array[:,0:N]
ds_array = ds_array[0:M,:]
im_geotrans = ds.GetGeoTransform()
cellsize = 1#im_geotrans[1]#栅格大小



def new4_3d(transmit,i,j,array):
    #发射点到接收点两点间水平距离
    D = math.sqrt(math.pow((i-transmit[0])*cellsize, 2)+math.pow((j-transmit[1])*cellsize, 2))
    #发射点到接收点两点间空间距离
    SD = math.sqrt(D*D+math.pow((array[i][j]+H1-array[transmit[0]][transmit[1]]-H), 2))
    #椭球短半轴长度
    HMi = 0.5*math.sqrt(SD*wavelength)
    #椭球长轴长度（将发射点和接收点作为椭球焦点）
    HMa = 0.5*math.sqrt(SD*SD+HMi*HMi)
    #投影最大椭圆长轴与半轴    HMi_2d = HMi

    HMa_2d = 0.5*math.sqrt(D*D+HMi_2d*HMi_2d)   

    sin = (array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
    cos = D/SD

    sin_xy = (j-transmit[1])*cellsize/D
    cos_xy = (i-transmit[0])*cellsize/D
    view = 0

    for d in range(0,91,degree):
  
        HMi_2d = HMi*math.cos(d*math.pi/180)

        if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
            s = abs(j-transmit[1])*cellsize/D #正弦值
            c = abs(i-transmit[0])*cellsize/D #余弦值
        else :
            s = abs(j-transmit[1])*cellsize/D #正弦值
            c = -abs(i-transmit[0])*cellsize/D #余弦值

        ls1 = list()
        ls2 = list()
        #先按行求
        if i < transmit[0]:
            h1 = i
            h2 = transmit[0]
        else :
            h1 = transmit[0]
            h2 = i
        #投影椭圆公式为：(x*c+y*s)**2/a**2+(x*s-y*c)**2/b**2=1
        for h in range(h1,h2):
            x = (h-((i+transmit[0])/2))*cellsize
            a = HMi_2d**2*s**2+HMa_2d**2*c**2
            b = (2*HMi_2d**2*c*s-2*HMa_2d**2*c*s)*x
            temp_c = (HMi_2d**2*c**2+HMa_2d**2*s**2)*x**2-HMa_2d**2*HMi_2d**2
            if b**2-4*a*temp_c >= 0:
                y1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                y2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a)
                y1 = (y1 + ((j+transmit[1])/2)*cellsize)/cellsize
                y2 = (y2 + ((j+transmit[1])/2)*cellsize)/cellsize
                #print(y1)
                #print(y2)
                if y1>=0 and y2>=0 and y1<M and y2<M:
                    if (transmit[1]-j)*(h-i)-(transmit[0]-i)*(y1-j)>0:
                        ls1.append(h)
                        ls1.append(y1)
                        ls2.append(h)
                        ls2.append(y2)
                    else:
                        ls2.append(h)
                        ls2.append(y1)
                        ls1.append(h)
                        ls1.append(y2)
                if y1<0 and y2>0 and y1<M and y2<M:
                    if (transmit[1]-j)*(h-i)-(transmit[0]-i)*(y2-j)>0:
                        ls1.append(h)
                        ls1.append(y2)
                    else :
                        ls2.append(h)
                        ls2.append(y2)
                if y2<0 and y1>0 and y1<M and y2<M:
                    if (transmit[1]-j)*(h-i)-(transmit[0]-i)*(y1-j)>0:
                        ls1.append(h)
                        ls1.append(y1)
                    else :
                        ls2.append(h)
                        ls2.append(y1)

        if j < transmit[1]:
            l1 = j
            l2 = transmit[1]
        else :
            l1 = transmit[1]
            l2 = j
        for l in range(l1,l2):
            y = (l-((j+transmit[1])/2))*cellsize
            a = HMi_2d**2*c**2+HMa_2d**2*s**2
            b = (2*HMi_2d**2*c*s-2*HMa_2d**2*c*s)*y
            temp_c = (HMi_2d**2*s**2+HMa_2d**2*c**2)*y**2-HMa_2d**2*HMi_2d**2
            if b**2-4*a*temp_c >= 0:
                x1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                x2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a)
                #print(x1)
                #print(x2)
                x1 = (x1 + ((i+transmit[0])/2)*cellsize)/cellsize
                x2 = (x2 + ((i+transmit[0])/2)*cellsize)/cellsize
 if x1>=0 and x2>=0 and x1 <N and x2 <N:
                    if (transmit[1]-j)*(x1-i)-(transmit[0]-i)*(l-j)>0:#（x1,l）在直线上方
                        ls1.append(x1)
                        ls1.append(l)
                        ls2.append(x2)
                        ls2.append(l)
                    else:
                        ls2.append(x1)
                        ls2.append(l)
                        ls1.append(x2)
                        ls1.append(l)
                if x1<0 and x2>0 and x1 <N and x2 <N:
                    if (transmit[1]-j)*(x2-i)-(transmit[0]-i)*(l-j)>0:
                        ls1.append(x2)
                        ls1.append(l)
                    else :
                        ls2.append(x2)
                        ls2.append(l)
                if x2<0 and x1>0 and x1 <N and x2 <N:
                    if (transmit[1]-j)*(x1-i)-(transmit[0]-i)*(l-j)>0:
                        ls1.append(x1)
                        ls1.append(l)
                    else :
                        ls2.append(x1)
                        ls2.append(l)
        print(ls1)
        print(ls2)

        if c>=0:
            if (array[i][j]+H1>=array[transmit[0]][transmit[1]]+H and i>=transmit[0]) \
                or (array[i][j]+H1<=array[transmit[0]][transmit[1]]+H and i<=transmit[0]):
                    s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                    c1 = abs(D)/SD
            else:
                s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                c1 = -abs(D)/SD
        else :
            if (array[i][j]+H1>=array[transmit[0]][transmit[1]]+H and i>=transmit[0]) \
                or (array[i][j]+H1<=array[transmit[0]][transmit[1]]+H and i<=transmit[0]):
                    s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                    c1 = -abs(D)/SD
            else:
                s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                c1 = abs(D)/SD
        
        flag = 0  

        if len(ls1)>0:
            #(ls1[k],ls1[k+1])代表当前所判断的单元
            for k in range(0,len(ls1),2):
 
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    if i<=transmit[0] and j<=transmit[1]:
                        temp_p = ls1[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls1[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                else:
                    if i<transmit[0] and j>transmit[1]:
                        temp_p = ls1[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls1[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                #再旋转
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                else:
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                
     
                a = 0.36*HMi**2*s1**2+HMa**2*c1**2
                b = 0.72*HMi**2*c1*s1*p-2*HMa**2*c1*s1*p
                temp_c = 0.36*HMi**2*c1**2*p**2+0.36*HMa**2*q**2+HMa**2*s1**2*p**2-0.36*HMi**2*HMa**2
                temp = b**2-4*a*temp_c
                #print(temp) 
                if b**2-4*a*temp_c>=0:
                    r1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                    r2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a) 
                    z1 = r1+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    z2 = r2+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    #print(z1)
                    #print(z1)
                    if z1>=z2:
                        #超过菲涅尔区的上界
                        if array[int(ls1[k])][int(ls1[k+1])]>z1:
                            break
                    else:
                        if array[int(ls1[k])][int(ls1[k+1])]>z2:
                            break
                if k>=(len(ls1)-2):
                    view = view+1
                    if d == 90:
                        flag = flag+1
          
        if len(ls2)>0:
            for k in range(0,len(ls2),2):
      
                '''
                temp_p = (ls2[k]-(i+transmit[0])/2)*cellsize
                temp_q = (ls2[k+1]-(j+transmit[1])/2)*cellsize
                p = temp_p*c+temp_q*s
                q = temp_p*s-temp_q*c
                '''

                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    if i<=transmit[0] and j<=transmit[1]:
                        temp_p = ls2[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls2[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                else:
                    if i<transmit[0] and j>transmit[1]:
                        temp_p = ls2[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls2[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                #再旋转
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                else:
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                
                #椭球在z轴上平移后的值z=z+array[transmit[0]][transmit[1]]
                a = 0.36*HMi**2*s1**2+HMa**2*c1**2
                b = 0.72*HMi**2*c1*s1*p-2*HMa**2*c1*s1*p
                temp_c = 0.36*HMi**2*c1**2*p**2+0.36*HMa**2*q**2+HMa**2*s1**2*p**2-0.36*HMi**2*HMa**2
                temp = b**2-4*a*temp_c
                #print(temp)
                if b**2-4*a*temp_c>=0:
                    r1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                    r2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a) 
                    z1 = r1+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    z2 = r2+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    #print(z1)
                    #print(z2)
                    if z1>=z2:
                        #超过菲涅尔区的上界
                        if array[int(ls2[k])][int(ls2[k+1])]>z1:
                            break
                    else:
                        if array[int(ls2[k])][int(ls2[k+1])]>z2:
                            break
                if k>=(len(ls2)-2):
                    view = view+1
                    if d == 90:
                        flag = flag+1
        if flag == 2:
            view = view-1
    return view/(2*90/degree+1)

#将公式应用到基于椭圆曲线的方法上，得出真实的信号强度值
def new6_3d(transmit,i,j,array):
    #发射点到接收点两点间水平距离
    D = math.sqrt(math.pow((i-transmit[0])*cellsize, 2)+math.pow((j-transmit[1])*cellsize, 2))
    #发射点到接收点两点间空间距离
    SD = math.sqrt(D*D+math.pow((array[i][j]+H1-array[transmit[0]][transmit[1]]-H), 2))
    #椭球短半轴长度
    HMi = 0.5*math.sqrt(SD*wavelength)
    #椭球长轴长度（将发射点和接收点作为椭球焦点）
    HMa = 0.5*math.sqrt(SD*SD+HMi*HMi)
    #投影最大椭圆长轴与半轴
    HMi_2d = HMi
    #HMi_2d = 0.5*math.sqrt(wavelength*D)
    HMa_2d = 0.5*math.sqrt(D*D+HMi_2d*HMi_2d)   
    sin = (array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
    cos = D/SD

    sin_xy = (j-transmit[1])*cellsize/D
    cos_xy = (i-transmit[0])*cellsize/D
    view = 0

    total_i = 0

    total_j = 0

    flag = 0

    flag1 = 0
 
    flag2 = 0

    for d in range(0,91,degree):

        HMi_2d = HMi*math.cos(d*math.pi/180)

        if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
            s = abs(j-transmit[1])*cellsize/D #正弦值
            c = abs(i-transmit[0])*cellsize/D #余弦值
        else :
            s = abs(j-transmit[1])*cellsize/D #正弦值
            c = -abs(i-transmit[0])*cellsize/D #余弦值


        ls1 = list()#存放发射点到接收端直线上方的椭圆边界点
        ls2 = list()#存放发射点到接收端直线下方的椭圆边界点
        #先按行求
        if i < transmit[0]:
            h1 = i
            h2 = transmit[0]
        else :
            h1 = transmit[0]
            h2 = i
        #投影椭圆公式为：(x*c+y*s)**2/a**2+(x*s-y*c)**2/b**2=1
        for h in range(h1,h2):
            x = (h-((i+transmit[0])/2))*cellsize
            a = HMi_2d**2*s**2+HMa_2d**2*c**2
            b = (2*HMi_2d**2*c*s-2*HMa_2d**2*c*s)*x
            temp_c = (HMi_2d**2*c**2+HMa_2d**2*s**2)*x**2-HMa_2d**2*HMi_2d**2
            if b**2-4*a*temp_c >= 0:
                y1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                y2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a)
                y1 = (y1 + ((j+transmit[1])/2)*cellsize)/cellsize
                y2 = (y2 + ((j+transmit[1])/2)*cellsize)/cellsize
                #print(y1)
                #print(y2)
                if y1>=0 and y2>=0 and y1<M and y2<M:
                    if (transmit[1]-j)*(h-i)-(transmit[0]-i)*(y1-j)>0:#（h,y1）在直线上方
                        ls1.append(h)
                        ls1.append(y1)
                        ls2.append(h)
                        ls2.append(y2)
                    else:
                        ls2.append(h)
                        ls2.append(y1)
                        ls1.append(h)
                        ls1.append(y2)
                if y1<0 and y2>0 and y1<M and y2<M:
                    if (transmit[1]-j)*(h-i)-(transmit[0]-i)*(y2-j)>0:
                        ls1.append(h)
                        ls1.append(y2)
                    else :
                        ls2.append(h)
                        ls2.append(y2)
                if y2<0 and y1>0 and y1<M and y2<M:
                    if (transmit[1]-j)*(h-i)-(transmit[0]-i)*(y1-j)>0:
                        ls1.append(h)
                        ls1.append(y1)
                    else :
                        ls2.append(h)
                        ls2.append(y1)
        #再按列求
        if j < transmit[1]:
            l1 = j
            l2 = transmit[1]
        else :
            l1 = transmit[1]
            l2 = j
        for l in range(l1,l2):
            y = (l-((j+transmit[1])/2))*cellsize
            a = HMi_2d**2*c**2+HMa_2d**2*s**2
            b = (2*HMi_2d**2*c*s-2*HMa_2d**2*c*s)*y
            temp_c = (HMi_2d**2*s**2+HMa_2d**2*c**2)*y**2-HMa_2d**2*HMi_2d**2
            if b**2-4*a*temp_c >= 0:
                x1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                x2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a)
                #print(x1)
                #print(x2)
                x1 = (x1 + ((i+transmit[0])/2)*cellsize)/cellsize
                x2 = (x2 + ((i+transmit[0])/2)*cellsize)/cellsize

                if x1>=0 and x2>=0 and x1 <N and x2 <N:
                    if (transmit[1]-j)*(x1-i)-(transmit[0]-i)*(l-j)>0:
                        ls1.append(x1)
                        ls1.append(l)
                        ls2.append(x2)
                        ls2.append(l)
                    else:
                        ls2.append(x1)
                        ls2.append(l)
                        ls1.append(x2)
                        ls1.append(l)
                if x1<0 and x2>0 and x1 <N and x2 <N:
                    if (transmit[1]-j)*(x2-i)-(transmit[0]-i)*(l-j)>0:
                        ls1.append(x2)
                        ls1.append(l)
                    else :
                        ls2.append(x2)
                        ls2.append(l)
                if x2<0 and x1>0 and x1 <N and x2 <N:
                    if (transmit[1]-j)*(x1-i)-(transmit[0]-i)*(l-j)>0:
                        ls1.append(x1)
                        ls1.append(l)
                    else :
                        ls2.append(x1)
                        ls2.append(l)

        if c>=0:
            if (array[i][j]+H1>=array[transmit[0]][transmit[1]]+H and i>=transmit[0]) \
                or (array[i][j]+H1<=array[transmit[0]][transmit[1]]+H and i<=transmit[0]):
                    s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                    c1 = abs(D)/SD
            else:
                s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                c1 = -abs(D)/SD
        else :
            if (array[i][j]+H1>=array[transmit[0]][transmit[1]]+H and i>=transmit[0]) \
                or (array[i][j]+H1<=array[transmit[0]][transmit[1]]+H and i<=transmit[0]):
                    s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                    c1 = -abs(D)/SD
            else:
                s1 = abs(array[i][j]+H1-array[transmit[0]][transmit[1]]-H)/SD
                c1 = abs(D)/SD

        flag = 0

        flag1 = 0

        flag2 = 0
        
        #(x*c+y*s)**2/a**2+(x*s-y*c)**2/b**2=1

        if len(ls1)>0:
            for k in range(0,len(ls1),2):
               
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    if i<=transmit[0] and j<=transmit[1]:
                        temp_p = ls1[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls1[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                else:
                    if i<transmit[0] and j>transmit[1]:
                        temp_p = ls1[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls1[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls1[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                #再旋转
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                else:
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                
                #椭球在z轴上平移后的值z=z+array[transmit[0]][transmit[1]]
                a = 0.36*HMi**2*s1**2+HMa**2*c1**2
                b = 0.72*HMi**2*c1*s1*p-2*HMa**2*c1*s1*p
                temp_c = 0.36*HMi**2*c1**2*p**2+0.36*HMa**2*q**2+HMa**2*s1**2*p**2-0.36*HMi**2*HMa**2
                temp = b**2-4*a*temp_c
                #print(temp) 
                if b**2-4*a*temp_c>=0:
                    r1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                    r2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a) 
                    z1 = r1+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    z2 = r2+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    #print(z1)
                    #print(z1)
                    if z1>=z2:
                        #超过菲涅尔区的上界
                        if array[int(ls1[k])][int(ls1[k+1])]>z2:
                            flag1 = 1
                        if array[int(ls1[k])][int(ls1[k+1])]>z1/0.6:
                            flag2 = 1
                        
                    else:
                        if array[int(ls1[k])][int(ls1[k+1])]>z1:
                            flag1 = 1
                        if array[int(ls1[k])][int(ls1[k+1])]>z2/0.6:
                            flag2 = 1
                        
            if flag1==1:
                total_i = total_i+1
            if flag2==1:
                total_j = total_j+1
            if d == 90 and flag1==1:
                flag = flag+1
                total_i = total_i-1
        #重置为0进行下一步判断
        #flag用于判断当前是否为视线
        flag = 0
        #flag1用于判断当前点是否挡住第一菲涅尔区60%
        flag1 = 0
        #flag2用于判断当前点是否挡住第一菲尼尔区的上界
        flag2 = 0  
        if len(ls2)>0:
            for k in range(0,len(ls2),2):
                #经过平移过后的坐标点（p,q）,可直接代入椭球进行计算，椭球在xy轴无需再平移
                '''
                temp_p = (ls2[k]-(i+transmit[0])/2)*cellsize
                temp_q = (ls2[k+1]-(j+transmit[1])/2)*cellsize
                p = temp_p*c+temp_q*s
                q = temp_p*s-temp_q*c
                '''
                #斜率为正,平移
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    if i<=transmit[0] and j<=transmit[1]:
                        temp_p = ls2[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls2[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                else:
                    if i<transmit[0] and j>transmit[1]:
                        temp_p = ls2[k]*cellsize-transmit[0]*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-transmit[1]*cellsize+(HMa_2d-D/2)*s
                    else:
                        temp_p = ls2[k]*cellsize-i*cellsize+(HMa_2d-D/2)*c
                        temp_q = ls2[k+1]*cellsize-j*cellsize+(HMa_2d-D/2)*s
                #再旋转
                if (i>=transmit[0] and j>=transmit[1]) or (i<=transmit[0] and j<=transmit[1]):
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                else:
                    p = temp_p*c+temp_q*s
                    p = p-HMa_2d
                    q = -temp_p*s+temp_q*c
                
                #椭球在z轴上平移后的值z=z+array[transmit[0]][transmit[1]]
                a = 0.36*HMi**2*s1**2+HMa**2*c1**2
                b = 0.72*HMi**2*c1*s1*p-2*HMa**2*c1*s1*p
                temp_c = 0.36*HMi**2*c1**2*p**2+0.36*HMa**2*q**2+HMa**2*s1**2*p**2-0.36*HMi**2*HMa**2
                temp = b**2-4*a*temp_c
                #print(temp)
                if b**2-4*a*temp_c>=0:
                    r1=(-b+math.sqrt(b**2-4*a*temp_c))/(2*a)
                    r2=(-b-math.sqrt(b**2-4*a*temp_c))/(2*a) 
                    z1 = r1+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    z2 = r2+(array[i][j]+H1+array[transmit[0]][transmit[1]]+H)/2
                    #print(z1)
                    #print(z2)
                    if z1>=z2:
                        #超过菲涅尔区的上界
                        if array[int(ls2[k])][int(ls2[k+1])]>z2:
                            flag1 = 1
                        if array[int(ls2[k])][int(ls2[k+1])]>z1/0.6:
                            flag2 = 1
                    else:
                        if array[int(ls2[k])][int(ls2[k+1])]>z1:
                            flag1 = 1
                        if array[int(ls2[k])][int(ls2[k+1])]>z2/0.6:
                            flag2 = 1
            if flag1==1:
                total_i = total_i+1
            if flag2==1:
                total_j = total_j+1
            if d == 90 and flag1==1:
                flag = flag+1
                total_i = total_i-1
                
    #发射强度为21.89
    db = 21.89
    #发射频率
    Hz = 2330
    n = (2*90/degree+1)
    if total_j>0:
        total_j = total_j-1
    #print(total_i)
    #print(total_j)
    #遮挡产生的衰减
    attenuation1 = 35*flag/n+20*total_i/n+10*total_j/n
    #距离产生的衰减
    attenuation2 = 20*math.log(Hz,10)+20*math.log(SD/1000,10)+32.4    
    return db-attenuation1-attenuation2

    
def main():

    #print(ds_array.shape)
    
    #print(band1_array.shape)
    #print(band1_array)

    
    print(im_geotrans)
    print(ds_array)
    #指定接收区域范围
    #ds_receive = ds_array[:,200:N]
    #ds_receive = ds_receive[0:200,:]
    #print(ds_receive)

    #选定一个信号发射点
    ds_transmit = [0,0]
    print(ds_array[ds_transmit[0]][ds_transmit[1]])
    result = [[-1 for j in range(0, M)] for i in range(0, N)]

    #循环计算
    for i in range(0,N):
        for j in range(0,M):
            if i == ds_transmit[0] and j == ds_transmit[1]:
                result[i][j] = 0
                continue
            '''
            if (i == 1 and j== 0) or (i==1 and j==1) or (i==0 and j==1) or(i == 2 and j == 2):
                result[i][j] = 3
                continue
            '''

            new4_value = new4_3d(ds_transmit,i,j,ds_array)
            new5_value = new5_3d(ds_transmit,i,j,ds_array)
            new6_value = new6_3d(ds_transmit,i,j,ds_array)
            result[i][j] = new4_value
    print(result)
    count = 0
    for q in range(0,N):
        for w in range(0,M):
            if result[q][w]==0:
                count=count+1
    print(count)
    #创建图片
    fig = plt.figure()
    sub = fig.add_subplot(111)
    im = sub.imshow(result,cmap='jet',origin='lower')
    cb = fig.colorbar(im,extend = 'both')
    '''
    for i in range(0,M):
        for j in range(0,N):
            if result[j][i]>0 and result[j][i]<1:
                sub.text(i,j,'{:.2f}'.format(result[j][i]),ha='center',va='center')
            else:
                sub.text(i,j,'{:d}'.format(result[j][i]),ha='center',va='center')
    '''
    plt.show()
    
if __name__ == '__main__':
      main()         
      
      
      

           
    
    
    
        
    
    
        
        
        
        