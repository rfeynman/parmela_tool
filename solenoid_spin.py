'''
Created on Jul 5, 2018

@author: wange
'''
import numpy as np

from scipy import interpolate
from scipy.integrate import quad
import matplotlib.pyplot as plt


def readfile():
    with open("../BFIELD.dat") as b:
        bfile=[line.split() for line in b]
    b.close()
    b_usef=bfile[55:]
    b_headname=[elem for elem in b_usef[0]][1:]
    b_usef[0]=b_headname
    #print(b_usef[0])
    col_b={b_usef[0][i]:[line[i] for line in b_usef[1:]] for i in range(len(b_usef[0]))}

    z_magnet=np.array([float(i) for i in col_b['Z(cm)']])
    b_magnet=np.array([float(i) for i in col_b['B(V/cm)']])
    magnet_bvsz=interpolate.interp1d(z_magnet,b_magnet)
    #print(z_magnet)
    #b_usef=np.array(b_usef)
    #b_col=b_usef.transpose()
    #col_b=np.asarray(col_b)
    z_ini=z_magnet[0]
    z_end=z_magnet[-1]
   
    
    with open("../TIMESTEPEMITTANCE.dat") as t:
        tfile=[line.split() for line in t]
    t.close()
    t_usef=tfile[84:]
    t_headname=[elem for elem in t_usef[0]][1:]
    t_usef[0]=t_headname
    col_t={t_usef[0][i]:[line[i] for line in t_usef[1:]] for i in range(len(t_usef[0]))}
    z_beam=np.array([float(i) for i in col_t['Z(cm)']])
    energy_beam=np.array([float(i) for i in col_t['<kE>(MeV)']])
    gamma_beam=(energy_beam+0.511)/0.511
    p_beam=energy_beam*1000*((gamma_beam+1)/(gamma_beam-1))**0.5  #momentum of the beam in unit keV/c
    beam_pvsz=interpolate.interp1d(z_beam,p_beam)
    #print(z_ini,z_end, magnet_bvsz(z_ini+3))
    #print('\n',z_beam,'\n',energy_beam)
    return z_magnet, b_magnet,z_beam,energy_beam,magnet_bvsz,beam_pvsz,z_ini,z_end
    
def spin(): #0.3(B*L)/p
    z_magnet, b_magnet,z_beam,energy_beam,magnet_bvsz,beam_pvsz,z_ini,z_end=readfile()
    print(z_ini,z_end,magnet_bvsz(z_ini+3),beam_pvsz(z_ini+3))
    #spin_phi=0.3*quad(lambda x:magnet_bvsz(x)/beam_pvsz(x),z_ini,z_end)
    nstep=200
    dstep=(z_end-z_ini)/nstep
    phase=[0]
    z=[z_ini]
    #print(phase[-1])
    #print(range(nstep))
    for i in range(nstep):
        stepphase=magnet_bvsz(z_ini+i*dstep)/beam_pvsz(z_ini+i*dstep)*dstep+phase[-1]
        phase.append(stepphase)
        z.append(z_ini+i*dstep)
    #print(phase)
    spin_phase=np.array(phase)
    spin_z=np.array(z)
    spin_phasevsz=np.column_stack((spin_z,spin_phase))
    print(spin_phasevsz) 
    return spin_phasevsz

def figure_plot():
    z_magnet, b_magnet,z_beam,energy_beam,magnet_bvsz,beam_pvsz,z_ini,z_end=readfile()
    spin_phasevsz=spin()
    fig,(ax1,ax3)=plt.subplots(2,1,sharex=True,figsize=(16,8))
    ax2=ax1.twinx()
    ax1.plot(z_magnet,b_magnet,'g-')
    ax2.plot(z_magnet,beam_pvsz(z_magnet),'b-')
    
    #plt.xlim(1.6,6.5)
    #ticks=np.arange(start-0.2,stop+0.2,0.2)
    #ax1.set_xticks(ticks)
    #print ticks
    ax1.set_ylabel('magnetic field[Gs]',color='g')
   
    ax2.set_ylabel('momentum[keV/c]',color='b')

    ax2.margins(0.02) 
    print(spin_phasevsz[:,0])
    ax3.set_xlabel('z[m]')
    ax3.set_ylabel('spin angle[deg]',color='r')
    ax3.plot(spin_phasevsz[:,0],spin_phasevsz[:,1],'r.')
    ax1.set_title('Electron beam momentum in the continuous focusing channel ')
    ax3.set_title('Spin phase angle envolution in focusing channel')
    plt.savefig("spin_angle_atfocusing.png")
    plt.show() 
    return 
    
if __name__ == '__main__':
    figure_plot()