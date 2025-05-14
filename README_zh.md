# kekule_lattice

## 1 安装教程

1.  下载库：pybinding, numpy

## 2 函数说明

创建一个 O 型 kekule 晶格

```
O_keku(onsite_energy=[POE,POE,POE,NOE,NOE,NOE])
```

裁切一个正六边形
其中 o 为原点坐标，r 为六边形外切圆半径

```
trapezoid(o=[0,0], r=10)
```

计算波函数

```
probability(model)
```

## 3 使用教程

### 3.1 设置程序

1. 打开 Calculation.py
2. 在代码的起始部分为变量设置，请先对其切割形状的边缘长度等物理量进行调整
3. 运行程序，观察程序给出的切割图片，若切割不合法，则停止程序重新设置参数，若无误，进入下一步
4. 等待程序运行，可以通过下方进度条来了解程序的预估运行时间
5. 进入 Calculation.py 所在的文件夹，打开 result 文件夹，里面存放了本次运行的全部结果
