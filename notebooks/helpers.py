import numpy as np
import matplotlib as mpl

rg_nodes = np.array([[4700, 3.7],
                     [3000, 0.],
                     [4250, 0.],
                     [5600, 2.3],
                     [5800, 3.3],
                     [4700, 3.7]])

ms_nodes = np.array([[3700, 5.2],
                     [3700, 4.5],
                     [4300, 4.2],
                     [5100, 4.2],
                     [8000, 3.4],
                     [8000, 4.75],
                     [3700, 5.2]])

sg_nodes = np.array([[4700, 3.7],
                     [5600, 3.4],
                     [6000, 3.8],
                     [5100, 4.2],
                     [4700, 3.7]])

ms_path = mpl.path.Path(ms_nodes)
rg_path = mpl.path.Path(rg_nodes)
sg_path = mpl.path.Path(sg_nodes)

def fast_m2_min(m1, mf):
    return (2*mf + np.power(2,0.6666666666666666)*
      np.power(mf*(27*pow(m1,2) + 18*m1*mf + 2*np.power(mf,2) + 
          3*np.sqrt(3)*np.sqrt(np.power(m1,3)*(27*m1 + 4*mf))),0.3333333333333333) + 
     (2*mf*(6*m1 + mf))/
      np.power((27*np.power(m1,2)*mf)/2. + 9*m1*np.power(mf,2) + np.power(mf,3) + 
        (3*np.sqrt(3)*mf*np.sqrt(np.power(m1,3)*(27*m1 + 4*mf)))/2.,0.3333333333333333))/6.