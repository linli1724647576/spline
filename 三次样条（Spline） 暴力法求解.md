## 三次样条（Spline） 暴力法求解

### 原理

<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308122145012.png" alt="image-20200308122145012" style="zoom:50%;" />![image-20200308122205949](C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308122205949.png)

<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308122145012.png" alt="image-20200308122145012" style="zoom:50%;" />![image-20200308122205949](C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308122205949.png)

<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308122243438.png" alt="image-20200308122243438" style="zoom:50%;" />



### 求解方法

本次采用暴力法求解，三次样条函数每一段是一个三次方程，设为 ax<sup>3</sup> + bx<sup>2</sup> +cx+d。  假设取11个待插值的点，这里就有十段函数。总共有40个未知数，需要40个方程求解。

首先输入方程两边相邻节点处函数值相等的方程为2n-2个方程，共有**18**个方程，然后输入端点处的函数值，有**2**个方程，端点函数一阶导数值相等为n-1个方程，这里有**9**个方程，端点函数二阶导数值相等，有n-1个方程，这里有**9**个方程，还差两个方程，需要人为补充条件。

情况一是端点处的函数值的二阶导数给出，为**2**个方程（默认为0）

情况二是端点处的函数值的一阶导数给出，为**2**个方程

所以总共的方程数就为**40**个，可以解出所有系数



### 算法解析

<table>
  <tr>
    <th width=40%, bgcolor=yellow>函数</th>
    <th width="60%", bgcolor=yellow>解释</th>
  </tr>
  <tr>
    <td> __init__(self,x,y,k)  </td>
    <td> 初始化函数，给定x,y,k，k为补充方程的情况  </td>
  </tr>
  <tr>
    <td> draw_pic(self,words, x, y) </td>
    <td> 绘图函数 </td>
  <tr>
    <td> spline3_Parameters(self,x_vec,k) </td>
    <td>  求解方程系数矩阵 </td>
  </tr>
  <tr>
    <td> solution_of_equation(self,parametes, x) </td>
    <td>  计算样条函数的系数 </td>
  </tr>
  <tr>
    <td> calculate(self,paremeters, x) </td>
    <td>  根据所给参数，计算三次函数的函数值 </td>
  </tr>
  <tr>
    <td> final(self) </td>
    <td>  输出图像 </td>
  </tr>    
</table>



### 源代码

```python
import math
import numpy as np
import matplotlib.pyplot  as plt
from sympy import *
from pylab import mpl

class Spline3(object):
    def __init__(self,x,y,k):
        self.x_init4 = x
        self.fucy = y
        self.k = k
        
    #画图
    def draw_pic(self,words, x, y):
        fig = plt.figure()
        plt.plot(x, y, label='interpolation')
        #plt.plot(x, self.fucy, label='raw')
        plt.scatter(self.x_init4,self.fucy,label='scatter',color='red')
        plt.legend()
        plt.title(words, FontProperties='SimHei')
        plt.show()
    #     plt.savefig(words + '.png')
    #     plt.close(fig)
    #     plt.show()
        pass
    
    def spline3_Parameters(self,x_vec,k):
        # parameter为二维数组，用来存放参数，size_of_Interval是用来存放区间的个数
        x_new = np.array(x_vec)
        parameter = []
        size_of_Interval = len(x_new) - 1;
        i = 1
        # 首先输入方程两边相邻节点处函数值相等的方程为2n-2个方程
        while i < len(x_new) - 1:
            data = np.zeros(size_of_Interval * 4)
            data[(i - 1) * 4] = x_new[i] * x_new[i] * x_new[i]
            data[(i - 1) * 4 + 1] = x_new[i] * x_new[i]
            data[(i - 1) * 4 + 2] = x_new[i]
            data[(i - 1) * 4 + 3] = 1

            data1 = np.zeros(size_of_Interval * 4)
            data1[i * 4] = x_new[i] * x_new[i] * x_new[i]
            data1[i * 4 + 1] = x_new[i] * x_new[i]
            data1[i * 4 + 2] = x_new[i]
            data1[i * 4 + 3] = 1

            parameter.append(data)
            parameter.append(data1)
            i += 1
        # 输入端点处的函数值。为两个方程, 加上前面的2n - 2个方程，一共2n个方程
        data = np.zeros(size_of_Interval * 4)
        data[0] = x_new[0] * x_new[0] * x_new[0]
        data[1] = x_new[0] * x_new[0]
        data[2] = x_new[0]
        data[3] = 1
        parameter.append(data)

        data = np.zeros(size_of_Interval * 4)
        data[(size_of_Interval - 1) * 4] = x_new[-1] * x_new[-1] * x_new[-1]
        data[(size_of_Interval - 1) * 4 + 1] = x_new[-1] * x_new[-1]
        data[(size_of_Interval - 1) * 4 + 2] = x_new[-1]
        data[(size_of_Interval - 1) * 4 + 3] = 1
        parameter.append(data)
        # 端点函数一阶导数值相等为n-1个方程。加上前面的方程为3n-1个方程。
        i = 1
        while i < size_of_Interval:
            data = np.zeros(size_of_Interval * 4)
            data[(i - 1) * 4] = 3 * x_new[i] * x_new[i]
            data[(i - 1) * 4 + 1] = 2 * x_new[i]
            data[(i - 1) * 4 + 2] = 1
            data[i * 4] = -3 * x_new[i] * x_new[i]
            data[i * 4 + 1] = -2 * x_new[i]
            data[i * 4 + 2] = -1
            parameter.append(data)
            i += 1
        # 端点函数二阶导数值相等为n-1个方程。加上前面的方程为4n-2个方程。
        i = 1
        while i < len(x_new) - 1:
            data = np.zeros(size_of_Interval * 4)
            data[(i - 1) * 4] = 6 * x_new[i]
            data[(i - 1) * 4 + 1] = 2
            data[i * 4] = -6 * x_new[i]
            data[i * 4 + 1] = -2
            parameter.append(data)
            i += 1
        #补充方程：
        if(k==1):
            # 端点处的函数值的二阶导数给出，为两个方程。总共为4n个方程。
            data = np.zeros(size_of_Interval * 4)
            data[0] = 6 * x_new[0]
            data[1] = 2
            parameter.append(data)
            data = np.zeros(size_of_Interval * 4)
            data[-4] = 6 * x_new[-1]
            data[-3] = 2
            parameter.append(data)
        if(k==2):
            #端点处的函数值的一阶导数给出，为两个方程。总共为4n个方程。
            data = np.zeros(size_of_Interval * 4)
            data[0] = 3 * x_new[i] * x_new[i]
            data[1] = 2 * x_new[i]
            data[2] = 1
            parameter.append(data)
            data = np.zeros(size_of_Interval * 4)
            data[-4] = 3 * x_new[i] * x_new[i]
            data[-3] = 2 * x_new[i]
            data[-2] = 1
            parameter.append(data)
        # df = pd.DataFrame(parameter)
        # df.to_csv('para.csv')
        #print(parameter)
        return parameter
    
    # 功能：计算样条函数的系数。
    # 参数：parametes为方程的系数，y为要插值函数的因变量。
    # 返回值：三次插值函数的系数。
    def solution_of_equation(self,parametes, x):
        size_of_Interval = len(x) - 1;
        result = np.zeros(size_of_Interval * 4)
        i = 1
        while i < size_of_Interval:
            result[(i - 1) * 2] = self.fucy[i]
            result[(i - 1) * 2 + 1] = self.fucy[i]
            i += 1
        result[(size_of_Interval - 1) * 2] = self.fucy[0]
        result[(size_of_Interval - 1) * 2 + 1] = self.fucy[-1]
        result[-2] = 5 / 13
        result[-1] = -5 / 13
        a = np.array(self.spline3_Parameters(x,self.k))
        b = np.array(result)
        # print(b)
    #     print(np.linalg.solve(a, b))
        return np.linalg.solve(a, b)


    # 功能：根据所给参数，计算三次函数的函数值：
    # 参数:parameters为二次函数的系数，x为自变量
    # 返回值：为函数的因变量
    def calculate(self,paremeters, x):
        result = []
        for data_x in x:
            result.append(
                paremeters[0] * data_x * data_x * data_x + paremeters[1] * data_x * data_x + paremeters[2] * data_x +
                paremeters[3])
        return result
    def final(self):
        result = self.solution_of_equation(self.spline3_Parameters(x_init4,self.k), self.x_init4)
        # print(spline3_Parameters(x_init4))
        # print(result)
        x_axis4 = []
        y_axis4 = []
        for i in range(10):
            temp = np.arange(-5 + i, -4 + i, 0.01)
            x_axis4 = np.append(x_axis4, temp)
            y_axis4 = np.append(y_axis4, self.calculate(
                [result[4 * i], result[1 + 4 * i], result[2 + 4 * i], result[3 + 4 * i]], temp))
        self.draw_pic('插值图', x_axis4, y_axis4)

        
x = np.arange(-5, 5.1, 1)
def func(y):
#     return 1 / (1 + y * y)
#     return np.cos(y)
#     return y
    return y**2
# x = np.arange(-5,5,1)
# y = [1,2,3,4,5,6,7,8,9,10,11]
#     return y**3
#     return np.sin(y)
y = func(x)
# y = np.array([1,2,3,4,5,6,7,8,9,12,5])
# print(y)
a = Spline3(x,y,1)
a.final()
```





### API

给定待插值节点x,和函数值y，参数k=1为一阶导数值给定，参数k=2为自然条件

示例

```python
x = np.arange(-5, 5.1, 1)
def func(y):
#     return 1 / (1 + y * y)
#     return np.cos(y)
    return y**2
y = func(x)
a = Spline3(x,y,1)
a.final()
```



### 效果图

![image-20200308124831760](C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308124831760.png)

