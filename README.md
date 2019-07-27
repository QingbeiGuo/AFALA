A Localization Method Avoiding Flip Ambiguities for micro-UAVs with Bounded Distance Measurement Errors
======
   Journal: IEEE Transactions on Mobile Computing

   @article{GuoZLKS18,
     title={A Localization Method Avoiding Flip Ambiguities for micro-UAVs with Bounded Distance Measurement Errors},
     author={Guo, Qingbei and Zhang, Yuan and Lloret, Jaime and Kantarci, Burak and Seah, Winston KG},
     journal={IEEE Transactions on Mobile Computing},
     volume={18},
     pages={1718--1730},
     year={2019},
     publisher={IEEE}
      doi={10.1109/TMC.2018.2865462},
    }

    Q. Guo, Y. Zhang, J. Lloret, B. Kantarci, and W. K. Seah, “A localization method avoiding flip ambiguities for micro-uavs with bounded distance measurement errors,” IEEE Transactions on Mobile Computing, vol. 18, pp. 1718–1730, 2019.

    论文 paper：https://arxiv.org/abs/1807.09590
    源码 code：Github: https://github.com/QingbeiGuo/AFALA.git.
    博客 blogs：https://blog.csdn.net/guoqwer/article/details/95322929

Abstract—Localization is a fundamental function in cooperative control of micro unmanned aerial vehicles (UAVs), but is easily affected by flip ambiguities because of measurement errors and flying motions. This study proposes a localization method that can avoid the occurrence of flip ambiguities in bounded distance measurement errors and constrained flying
motions; to demonstrate its efficacy, the method is implemented on bilateration and trilateration. For bilateration, an improved bi-boundary model based on the unit disk graph model is created to compensate for the shortage of distance constraints, and two boundaries are estimated as the communication range constraint. The characteristic of the intersections of the communication range and distance constraints is studied to present a unique localization criterion which can avoid the occurrence of flip ambiguities. Similarly, for trilateration, another unique localization criterion for avoiding flip ambiguities is proposed according to the characteristic of the intersections of three distance constraints. The theoretical proof shows that these proposed criteria are correct. A localization algorithm is constructed based on these two criteria. The algorithm is validated using simulations for different scenarios and parameters, and the proposed method is shown to provide excellent localization performance in terms of average estimated error. 

摘要：在无线传感器网络中定位是一个基本的功能，但是由于测量误差它容易受到翻转的影响。本研究提出一种能够在有限距离测量误差下避免翻转的定位方法，这一方法在两边定位和三边定位中进行了实现。对于两边定位，我们创建了一个基于单元圆盘图（unit disk graph）模型的双边界模型以弥补距离约束的不足，并且模型中的双边界被评估作为通讯范围的约束。我们研究了通讯范围和距离约束的交集特征提出了一个避免翻转的唯一性定位准则。类似的，对于三边定位，根据三边约束形成的交集特征提出了另一个避免翻转的唯一性定位准则。理论证明，我们提出的定位准则是正确的。我们基于这两个准则建立了一个定位算法，而且在不同的场景和参数下进行的实验评估，实验结果表明我们提出的算法能够完全避免翻转的发生，有效降低了节点的定位误差。


