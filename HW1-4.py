import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#=========================================================================================================================
#========================================Preparation======================================================================
#=========================================================================================================================

def u0(Xarray):
    u=np.array([])
    for x in Xarray:
        if ((x<=0)or (x>=2)):
            u=np.append(u,[0])
        if ((x>0) and (x<1)):
            u=np.append(u,[0.5*x])
        if ((x>1) and (x<2)):
            u=np.append(u,[1-0.5*x])

    return u

def u1(Xarray):  #I want to test with a smooth fonction
    return np.sin(1*Xarray)

def u2(Xarray):

    u=np.array([])
    for x in Xarray:
        if x<=1:
            u=np.append(u,[1])
        if x >1:
            u=np.append(u,[0])

    return u
        
            

def up_wind(U0,Mu,N,createMap=False):
    'If createMap==False, give last iteration, else, give whole map'
    i=0   #number of iteration
    Xnew=np.array(U0) #Copie of the array, not the array
    if (createMap):
        Map=np.array(U0)
    
    while (i<N):
        Xold=np.array(Xnew)
        j=1    #calculate j term of Xnew
        while (j<len(U0)):
            Xnew[j]=Xold[j]-Mu*(Xold[j]-Xold[j-1])
            j+=1
            
        Xnew[0]=Xnew[1]   #We cannot find this value, it is define by default here, that why we made a value change of Xa
        
        if (createMap):
            Map=np.append(Map,Xnew)      
        i+=1
    if (createMap):
        Map=np.reshape(Map,(N+1,len(U0)))
        return Map
    else:
        return Xnew

def mc_cormack(U0,Mu,N,createMap=False):
    Xnew=np.array(U0)
    if (createMap):
        Map=np.array(U0)
    i=0         #Number of iteration
    while i<N:
        Xold=np.array(Xnew)
        j=1         #Let's calculate term j of Xnew
        while (j<len(U0)-1):
            #First, intermediate value at (n+1/2): Xmiddle=[u^(n+1/2)_(j-1) ; u^(n+1/2)_(j)]
            Xmiddle=np.array([(Xold[j-1]-Mu*(Xold[j]-Xold[j-1])) , (Xold[j]-Mu*(Xold[j+1]-Xold[j]))  ])  
            #Then, Calculate Xnew[j]
            Xnew[j]=0.5*(Xold[j]+Xmiddle[1])-0.5*Mu*(Xmiddle[1]-Xmiddle[0])
            j+=1
        
        #We need To add Xnew[0]and Xnew[end]
        Xnew[0]=Xnew[1]
        Xnew[len(Xnew)-1:]=np.array([Xnew[len(Xnew)-2:-1]])
        
        if (createMap):
            Map=np.append(Map,Xnew)      
        i+=1     
               
    if (createMap):
        Map=np.reshape(Map,(N+1,len(U0)))
        return Map
    else:
        return Xnew

def create_array2D(Mu,N,Xa,Xb,Dx,Fonction):
    'Create an array with all the soluttion: Vector[a][k]-> a=0-> UpWind: a=1->MacCormack  k is which case (different mu) we use \
    fonction is which initial condition you what to '
    
    UPW=np.array([])                                                    #CleanUp
    MCC=np.array([])
    
    i=0
    while(i<len(Mu)):
        XaUpW=np.min(Xa-N[i]*Dx-Dx)  #we need to have point before the zone of interest
        XbUpW=Xb                   

        XaMcC=np.min(Xa-N[i]*Dx-Dx)  #we need to have point before the zone of interest
        XbMcC=np.max(Xb+N[i]*Dx+Dx)         
        
        X=np.arange(Xa,Xb,Dx)
        XUpW = np.arange(XaUpW, XbUpW, Dx)
        XMcC = np.arange(XaMcC, XbMcC, Dx)
        
        if(Fonction==0):
            U0UpW=u0(XUpW)
            U0McC=u0(XMcC)       
        
        if(Fonction==1):
            U0UpW=u1(XUpW)
            U0McC=u1(XMcC)
        
        if(Fonction==2):
            U0UpW=u2(XUpW)
            U0McC=u2(XMcC)
        
        
        solUpW=up_wind(U0UpW,Mu[i],N[i],False)
        solMcC=mc_cormack(U0McC,Mu[i],N[i],False)
        
        solUpW_Shaped = solUpW[N[i]+1:]                                 #reshaping ignoring extra X
        solMcC_Shaped=solMcC[N[i]+1:-N[i]-1]


        
        UPW=np.append(UPW,solUpW_Shaped)
        MCC=np.append(MCC,solMcC_Shaped)
        
        i+=1
    
    UPW=np.reshape(UPW,(len(Mu),len(X)))
    MCC=np.reshape(MCC,(len(Mu),len(X)))

    Vector=(UPW,MCC)
    return(Vector)

def final_plot(UPW_MCC,Xa,Xb,Dx,C,Mu,Tfinale,Fonction):
    X=np.arange(Xa,Xb,Dx)

    Name=np.array(['Exact solution','Mu=' + str(Mu[0]),'Mu=' + str(Mu[1]),'Mu=' + str(Mu[2]) ])

 
    plt.figure('Mac_Cormack')

    if Fonction==1:
        plt.plot(X,u1(X-C*Tfinale),'b',X,UPW_MCC[1][0],'-.',X,UPW_MCC[1][1],'-.',X,UPW_MCC[1][2],'-.')
    if Fonction==0:
        plt.plot(X,u0(X-C*Tfinale),'b',X,UPW_MCC[1][0],'-.',X,UPW_MCC[1][1],'-.',X,UPW_MCC[1][2],'-.')
    plt.legend( (Name[0],Name[1],Name[2],Name[3]) )
    plt.xlabel('x')
    plt.ylabel('Value at ' +str(Tfinale)+' s')
    plt.savefig('MCC.png')
    plt.show()
    
    plt.figure('Up Wind')

    if Fonction==1:
        plt.plot(X,u1(X-C*Tfinale),X,UPW_MCC[0][0],'-.',X,UPW_MCC[0][1],'-.',X,UPW_MCC[0][2],'-.')
    if Fonction==0:
        plt.plot(X,u0(X-C*Tfinale),'b',X,UPW_MCC[0][0],'-.',X,UPW_MCC[0][1],'-.',X,UPW_MCC[0][2],'-.')
    plt.legend( (Name[0],Name[1],Name[2],Name[3]) )
    plt.xlabel('x')
    plt.ylabel('Value at ' +str(Tfinale)+' s')
    plt.savefig('UPW.png')
    plt.show()

#=================================FINAL MAIN============================
def main():
    #Constant
    Tfinale=2 #final time
    Dx=0.2
    Mu=np.array([0.25,0.5,1.])
    C=0.5
    Dt=Dx*Mu/C
    N=np.int64(Tfinale/Dt) #number of iteration
    Xa,Xb= -2, 5 #Define our zone of interrest
    fonction=0
    
    UPW_MCC=create_array2D(Mu,N,Xa,Xb,Dx,fonction)
    final_plot(UPW_MCC,Xa,Xb,Dx,C,Mu,Tfinale,fonction)
    #~ print(np.argmax(UPW_MCC[0][0])*Dx-2,np.argmax(UPW_MCC[0][1]*Dx-2),np.argmax(UPW_MCC[1][0])*Dx-2)
    i=0

#=======================================================================

if __name__=='__main__':
    main()






#~ #================TEST MAIN==============================================
#~ def main():
    #~ #Constant
    #~ Tfinale=2 #final time
    #~ Dx=0.2
    #~ Mu=np.array([0.25,0.5,1.])
    #~ C=0.5
    #~ Dt=Dx*Mu/C
    #~ N=np.int64(Tfinale/Dt) #number of iteration
    #~ Xa,Xb= -2, 5 #Define our zone of interrest
    #~ fonction=0
    #~ X=np.arange(Xa,Xb,Dx)
    #~ U0=u0(X)
    #~ plt.plot(X,U0,X,up_wind(U0,1,1,False))
    #~ plt.show()
#~ #======================================================================= 

#=================================PLOT 3D===============================
#~ solUpW=sol[:,int(XaUpW-Xa):]
#~ solMcC=solMcC[:,N[0]+1:-N[0]-1]
#~ fig = plt.figure()
#~ ax = fig.add_subplot(111, projection='3d')


#~ x=XUpW = np.arange(Xa, Xb, Dx)
#~ t=np.arange(0,Tfinale+Dt[0],Dt[0])
#~ X,T=np.meshgrid(x,t)
#~ print(X.shape)
#~ print(solMcC.shape)

#~ ax.plot_wireframe(X, T, solMcC)

#~ plt.show()
#=======================================================================

    #~ plt.figure('Up_wind-0')
    #~ if Fonction==1:
        #~ plt.plot(X,u1(X-C*Tfinale),'b',X,UPW_MCC[0][0])
    #~ if Fonction==0:
        #~ plt.plot(X,u0(X-C*Tfinale),'b',X,UPW_MCC[0][0])
    #~ plt.savefig('UPW-0.png')
    #~ plt.legend( (Name[0],Name[1] ))
    #~ plt.xlabel('x')
    #~ plt.ylabel('Value at ' +str(Tfinale)+' s')
    #~ plt.show()
    
    #~ plt.figure('Up_wind-1')
    #~ if Fonction==1:
        #~ plt.plot(X,u1(X-C*Tfinale),'b',X,UPW_MCC[0][1])
    #~ if Fonction==0:
        #~ plt.plot(X,u0(X-C*Tfinale),'b',X,UPW_MCC[0][1])
    #~ plt.savefig('UPW-1.png')
    #~ plt.legend( (Name[0],Name[2]) )
    #~ plt.xlabel('x')
    #~ plt.ylabel('Value at ' +str(Tfinale)+' s')
    #~ plt.show()
    
    #~ plt.figure('Up_wind-2')
    #~ if Fonction==1:
        #~ plt.plot(X,u1(X-C*Tfinale),'b',X,UPW_MCC[0][2])
    #~ if Fonction==0:
        #~ plt.plot(X,u0(X-C*Tfinale),'b',X,UPW_MCC[0][2])
    #~ plt.savefig('UPW-2.png')
    #~ plt.legend( (Name[0],Name[3] ))
    #~ plt.xlabel('x')
    #~ plt.ylabel('Value at ' +str(Tfinale)+' s')
    #~ plt.show()
    
