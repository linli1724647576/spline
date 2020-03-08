

## 三次样条（Spline）三对角矩阵法求解

### 算法解析



样条插值公式：

<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308150726694.png" alt="image-20200308150726694" style="zoom: 67%;" />



已经y<sub>i</sub> ,y<sub>i+1</sub> ,  .......   m<sub>i</sub> ,m<sub>i+1</sub> 为代求系数...

<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308151034223.png" alt="image-20200308151034223" style="zoom:67%;" />



h<sub>i</sub>,h<sub>i+1</sub>,H<sub>i</sub>,H<sub>i+1</sub> 容易根据给出的点求得，主要要求 m<sub>i</sub> ,m<sub>i+1</sub> 

<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308151435926.png" alt="image-20200308151435926" style="zoom:50%;" />



<img src="C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308151513361.png" alt="image-20200308151513361" style="zoom:50%;" />



a. 自由边界(Natural)

![image-20200308151623649](C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308151623649.png)



b. 固定边界(Clamped)

![image-20200308151706120](C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308151706120.png)



两种条件下的mi的求解，参考博文https://blog.csdn.net/flyingleo1981/article/details/53008931



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
    <td> def spline1(self) </td>
    <td> 固定条件下的方程 </td>
  <tr>
    <td> def spline(self) </td>
    <td>  自然条件下的方程 </td>
  </tr>
  <tr>
    <td> def cal(self,xi, xii, i) </td>
    <td>  计算一段的函数 </td>
  </tr>
  <tr>
    <td> def calnf(self) </td>
    <td>  计算每一段的函数 </td>
  </tr>
  <tr>
    <td> def nfSub(self,x, nf) </td>
    <td>  求值 </td>
  </tr>    
  <tr>
    <td> def calf(self,f,x) </td>
    <td>  求值 </td>
  </tr>   
  <tr>
    <td> def draw(self,nf) </td>
    <td>  画图 </td>
  </tr>   
  <tr>
    <td> def final(self) </td>
    <td> 显示结果 </td>
  </tr>  
</table>



### 源代码

```python
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

class Spline(object):
    def __init__(self,x,y,k):
        self.x = x
        self.y= y
        self.k = k
        self.n = symbols('n')
    #求解三对角方程，固定边界，一阶导数已知
    def spline1(self):
    #取-5-5十一个点
        Mat = np.eye(11, 11) * 4
        ds = []
        for i in range(11):
            h = 1
            alaph = 0.5
            if i==0:
                beta = 0
            if 1<=i<10:
                beta = 6*(y[i+1]-2*y[i]+y[i-1])
            if i==10:
                beta = 0
            ds.append(beta)
            if i == 0:
                Mat[i][0] = 2
                Mat[i][1] = 0
            elif i == 10:
                Mat[i][9] = h
                Mat[i][10] = 2
            else:
                Mat[i][i - 1] = h
                Mat[i][i + 1] = alaph
        ds = np.mat(ds)
        Mat = np.mat(Mat)
        Ms = ds * Mat.I
        self.Ms = Ms.tolist()[0]
    #     return Ms.tolist()[0]

    #求解三对角方程，自然边界，二阶导数为0
    def spline(self):
        #取-5-5十一个点
    #     global x,y
        Mat = np.eye(11, 11) * 2
        ds = []
        for i in range(11):
            h = 1
            alaph = 0.5
            if i==0:
                beta = 0
            if 1<=i<10:
                beta = 3*(y[i+1]-y[i-1])/(2*h)
            if i==10:
                beta = 0
            ds.append(beta)
            if i == 0:
                Mat[i][0] = 1
                Mat[i][1] = alaph
            elif i == 10:
                Mat[i][9] = 1 - alaph
            else:
                Mat[i][i - 1] = 1 - alaph
                Mat[i][i]
                Mat[i][i + 1] = alaph
        ds = np.mat(ds)
        Mat = np.mat(Mat)
        Ms = ds * Mat.I
        self.Ms = Ms.tolist()[0]
#         return Ms.tolist()[0]
    #计算每一段的插值函数
    def cal(self,xi, xii, i):
#         Ms = self.spline()
        yi = self.y[i]
        yii = self.y[i+1]
        hi = (1+2*(self.n-xi)/(xii-xi))*((self.n-xii)/(xi-xii))**2
        hii = (1+2*(self.n-xii)/(xi-xii))*((self.n-xi)/(xii-xi))**2
        Hi = (self.n-xi)*((self.n-xii)/(xi-xii))**2
        Hii = (self.n-xii)*((self.n-xi)/(xii-xi))**2
        I = hi*yi+hii*yii+Hi*self.Ms[i]+Hii*self.Ms[i+1]
        return I
    def calnf(self):
        nf = []
        for i in range(len(self.x) - 1):
            nf.append(self.cal(self.x[i], self.x[i + 1], i))
        return nf
    #求值画图
    def nfSub(self,x, nf):
        tempx = np.array(range(11)) - 5
        dx = []
        for i in range(10):
            labelx = []
            for j in range(len(x)):
                if x[j] >= tempx[i] and x[j] < tempx[i + 1]:
                    labelx.append(x[j])
                elif i == 9 and x[j] >= tempx[i] and x[j] <= tempx[i + 1]:
                    labelx.append(x[j])
            dx = dx + self.calf(nf[i], labelx)
        return np.array(dx)
    def calf(self,f,x):
        y = []
        for i in x:
            y.append(f.subs(self.n, i))
        return y 
    #画图
    def draw(self,nf):
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        x = np.linspace(-5, 5, 101)
        Ly = self.nfSub(x,nf)
        plt.plot(x, Ly, label='三次样条插值函数')
        plt.scatter(self.x,self.y,label='scatter',color='red')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()

        plt.savefig('1.png')
        plt.show() 

    def final(self):
        init_printing(use_unicode=True)
        if(self.k==1):
            Ms = self.spline()
        elif(self.k==2):
            Ms = self.spline1()
        self.nf = self.calnf()
        self.draw(self.nf)
        

        
x = np.arange(-5, 5.1, 1)
def func(y):
#     return 1 / (1 + y * y)
    return np.cos(y)
#     return y
#     return y**3
y = func(x)

a = Spline(x,y,1)
a.final()
```





### API

给定待插值节点x,和函数值y，参数k=1为自然条件，参数k=2为固定条件

示例

```python
x = np.arange(-5, 5.1, 1)
def func(y):
#     return 1 / (1 + y * y)
#     return np.cos(y)
    return y**2
y = func(x)
a = Spline(x,y,1)
a.final()
```



### 效果图

![image-20200308124831760](C:\Users\HUAWEI\AppData\Roaming\Typora\typora-user-images\image-20200308124831760.png)

