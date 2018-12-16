# Change: Modifying so that the ends are straight coarse bricks
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
import csv
import os

def NodeGen2DV45(x0,xl,y0,yl,z0,elemLenX,elemLenY,numElemX,numElemY,shiftX,shiftY):
    # Nodal coordinated
    nodeX1=np.linspace(x0,xl,numElemX);
    nodeY1=y0+np.zeros(np.shape(nodeX1)) #np.arange(0,specLenY+elemLenY,elemLenY);
    #
    nodeX2=np.linspace(x0+shiftX,xl-shiftX,numElemX-1);
    nodeY2=y0+np.zeros(np.shape(nodeX2))+shiftY
    #
    # Create all nodes
    count=1;
    
    Node=np.array([[0,0,0,0]])
    for j in range(0,int(numElemY)-1):
        for i in range(0,len(nodeX1)):
            Node=np.append(Node,[[int(count+i),nodeX1[i],nodeY1[i]+j*elemLenY,z0]],axis=0)
        count=len(Node)
        
        for i in range(0,len(nodeX2)):
            Node=np.append(Node,[[int(count+i),nodeX2[i],nodeY2[i]+j*elemLenY,z0]],axis=0)
        count=len(Node)      
        
    # last line
    for i in range(0,len(nodeX1)):
        Node=np.append(Node,[[int(count+i),nodeX1[i],nodeY1[i]+(j+1)*elemLenY,z0]],axis=0)
            
    Node=Node[1:len(Node)]
    return Node

def NodeGen2DV90(x0,xl,y0,yl,z0,elemLenX,elemLenY,numElemX,numElemY):
    # Nodal coordinated
    nodeX1=np.linspace(x0,xl,numElemX);
    nodeY1=y0+np.zeros(np.shape(nodeX1)) #np.arange(0,specLenY+elemLenY,elemLenY)
    # Create all nodes
    count=1;
    
    Node=np.array([[0,0,0,0]])
    for j in range(0,int(numElemY)):
        for i in range(0,len(nodeX1)):
            Node=np.append(Node,[[int(count+i),nodeX1[i],nodeY1[i]+j*elemLenY,z0]],axis=0)
        count=len(Node)
#            
    Node=Node[1:len(Node)]
    elemLenX=nodeX1[1]-nodeX1[0]
    
    return Node,elemLenX,elemLenY

def FindNodes(loc,Node):
    NCorners=[[0,0,0,0]]
    for i in range(len(loc)):
        NCornersTmp=Node[(Node[:,1]==loc[i,0])]
        NCornersTmp=NCornersTmp[(NCornersTmp[:,2]==loc[i,1])]
        NCorners=np.append(NCorners,NCornersTmp, axis=0)
    NCorners=NCorners[1:len(NCorners)]
    return NCorners

def FindNodeRange(loc,Node,elemLenX,elemLenY):
    loc=[loc[0]-1e-5,loc[1]-1e-5]
    NCornersTmp=Node
    NCornersTmp=Node[(Node[:,1]>=loc[0])]
    NCornersTmp=NCornersTmp[(NCornersTmp[:,1]<=loc[0]+1.5*elemLenX)]
    NCornersTmp=NCornersTmp[(NCornersTmp[:,2]>=loc[1])]
    NCornersTmp=NCornersTmp[(NCornersTmp[:,2]<=loc[1]+1.5*elemLenY)]
    return NCornersTmp

def FindBoundariesV2(Node,x0,xl,y0,yl,numElemX,numElemY):
    # Find corners
    loc=np.array([[x0,y0],[xl,0],[0,yl],[xl,yl]])
    NCorners= FindNodes(loc,Node)
    
    # Find bottom edge
    Xrange=np.linspace(x0,xl,numElemX)
    Yrange=np.ones(np.shape(Xrange))*y0
    loc=np.transpose(np.array([Xrange,Yrange]))
    NBtmEdge= FindNodes(loc,Node)
    
    # Find top edge
    Xrange=np.linspace(x0,xl,numElemX)
    Yrange=np.ones(np.shape(Xrange))*yl
    loc=np.transpose(np.array([Xrange,Yrange]))
    NTopEdge= FindNodes(loc,Node)
    
    # Find left edge
    Yrange=np.linspace(y0,yl,numElemY)
    Xrange=np.ones(np.shape(Yrange))*x0
    loc=np.transpose(np.array([Xrange,Yrange]))
    NLeftEdge= FindNodes(loc,Node)
    
    # Find right edge
    Yrange=np.linspace(y0,yl,numElemY)
    Xrange=np.ones(np.shape(Yrange))*xl
    loc=np.transpose(np.array([Xrange,Yrange]))
    NRightEdge= FindNodes(loc,Node)
    
    NBoundary=np.append(NBtmEdge,NRightEdge,axis=0)
    NBoundary=np.append(NBoundary,NTopEdge,axis=0)
    NBoundary=np.append(NBoundary,NLeftEdge,axis=0)
    
    return NCorners,NBtmEdge,NTopEdge,NLeftEdge,NRightEdge,NBoundary

def FindBoundaries(Node,specLenX,specLenY,elemLenX,elemLenY):
    # Find corners
    loc=np.array([[0,0],[specLenX,0],[0,specLenY],[specLenX,specLenY]])
    NCorners= FindNodes(loc,Node)
    
    # Find bottom edge
    Xrange=np.arange(0,specLenX,elemLenX)
    Yrange=np.ones(np.shape(Xrange))*0
    loc=np.transpose(np.array([Xrange,Yrange]))
    NBtmEdge= FindNodes(loc,Node)
    
    # Find top edge
    Xrange=np.arange(0,specLenX,elemLenX)
    Yrange=np.ones(np.shape(Xrange))*specLenY
    loc=np.transpose(np.array([Xrange,Yrange]))
    NTopEdge= FindNodes(loc,Node)
    
    # Find left edge
    Yrange=np.arange(0,specLenY,elemLenY)
    Xrange=np.ones(np.shape(Yrange))*0
    loc=np.transpose(np.array([Xrange,Yrange]))
    NLeftEdge= FindNodes(loc,Node)
    
    # Find right edge
    Yrange=np.arange(0,specLenY,elemLenY)
    Xrange=np.ones(np.shape(Yrange))*specLenX
    loc=np.transpose(np.array([Xrange,Yrange]))
    NRightEdge= FindNodes(loc,Node)
    
    NBoundary=np.append(NBtmEdge,NRightEdge,axis=0)
    NBoundary=np.append(NBoundary,NTopEdge,axis=0)
    NBoundary=np.append(NBoundary,NLeftEdge,axis=0)
    
    return NCorners,NBtmEdge,NTopEdge,NLeftEdge,NRightEdge,NBoundary

def DefineElem2D45(Node,NBtmEdge,NTopEdge,NLeftEdge,NRightEdge,shiftX,shiftY):
    A=spatial.cKDTree(Node[:,1:3])
    
    # Find nearest
    XYPnt=np.array([0.0,0.0])
    ElemQuad=np.array([[0,0,0,0,0]])
    ElemPyrd=np.array([[0,0,0,0]])
    eleCount=1
    
    for i in range(0,len(NBtmEdge)):
        idx=np.ones([1,3])*-1
        XYPnt=NBtmEdge[i,1:3]
        distance1,idx1 = A.query([XYPnt[0]+shiftX,XYPnt[1]+shiftY],k=1,distance_upper_bound=2)
        distance2,idx2 = A.query([XYPnt[0],XYPnt[1]],k=1,distance_upper_bound=2)    
        distance3,idx3 = A.query([XYPnt[0]+2*shiftX,XYPnt[1]],k=1,distance_upper_bound=2)    
        idx=[idx1,idx2,idx3]
        idxTmp=np.unique(idx)
        if len(idxTmp)==3:
            ElemPyrd=np.append(ElemPyrd,[[eleCount,Node[idx[0],0],Node[idx[1],0],Node[idx[2],0] ]],axis=0)        
            eleCount=eleCount+1    
    
    for i in range(0,len(NTopEdge)):
        idx=np.ones([1,3])*-1
        XYPnt=NTopEdge[i,1:3]
        distance1,idx1 = A.query([XYPnt[0]+shiftX,XYPnt[1]-shiftY],k=1,distance_upper_bound=2)
        distance2,idx2 = A.query([XYPnt[0]+2*shiftX,XYPnt[1]],k=1,distance_upper_bound=2)
        distance3,idx3 = A.query([XYPnt[0],XYPnt[1]],k=1,distance_upper_bound=2)            
        idx=[idx1,idx2,idx3]
        idxTmp=np.unique(idx)
        if len(idxTmp)==3:
            ElemPyrd=np.append(ElemPyrd,[[eleCount,Node[idx[0],0],Node[idx[1],0],Node[idx[2],0] ]],axis=0)        
            eleCount=eleCount+1    
    
    for i in range(0,len(NLeftEdge)):
        idx=np.ones([1,3])*-1
        XYPnt=NLeftEdge[i,1:3]
        distance1,idx1 = A.query([XYPnt[0],XYPnt[1]],k=1,distance_upper_bound=2)        
        distance2,idx2 = A.query([XYPnt[0]+shiftX,XYPnt[1]+shiftY],k=1,distance_upper_bound=2)
        distance3,idx3 = A.query([XYPnt[0],XYPnt[1]+2*shiftY],k=1,distance_upper_bound=2)    
        idx=[idx1,idx2,idx3]
        idxTmp=np.unique(idx)
        if len(idxTmp)==3:
            ElemPyrd=np.append(ElemPyrd,[[eleCount,Node[idx[0],0],Node[idx[1],0],Node[idx[2],0] ]],axis=0)        
            eleCount=eleCount+1    
    
    for i in range(0,len(NRightEdge)):
        idx=np.ones([1,3])*-1
        XYPnt=NRightEdge[i,1:3]
        distance1,idx1 = A.query([XYPnt[0],XYPnt[1]+2*shiftY],k=1,distance_upper_bound=2)    
        distance2,idx2 = A.query([XYPnt[0]-shiftX,XYPnt[1]+shiftY],k=1,distance_upper_bound=2)
        distance3,idx3 = A.query([XYPnt[0],XYPnt[1]],k=1,distance_upper_bound=2)    
        idx=[idx1,idx2,idx3]
        idxTmp=np.unique(idx)
        if len(idxTmp)==3:
            ElemPyrd=np.append(ElemPyrd,[[eleCount,Node[idx[0],0],Node[idx[1],0],Node[idx[2],0] ]],axis=0)        
            eleCount=eleCount+1 
    
    for i in range(0,len(Node)):
        idx=np.ones([1,4])*-1
        XYPnt=Node[i,1:3]
        distance1,idx1 = A.query([XYPnt[0]+shiftX,XYPnt[1]-shiftY],k=1,distance_upper_bound=2)
        distance2,idx2 = A.query([XYPnt[0]+2*shiftX,XYPnt[1]],k=1,distance_upper_bound=2)
        distance3,idx3 = A.query([XYPnt[0]+shiftX,XYPnt[1]+shiftY],k=1,distance_upper_bound=2)
        distance4,idx4 = A.query([XYPnt[0],XYPnt[1]],k=1,distance_upper_bound=2)
        idx=[idx1,idx2,idx3,idx4]
        idxTmp=np.unique(idx)
        if len(idxTmp)==4:
            ElemQuad=np.append(ElemQuad,[[eleCount,Node[idx[0],0],Node[idx[1],0],Node[idx[2],0],Node[idx[3],0] ]],axis=0)
            eleCount=eleCount+1        
    
    ElemQuad=ElemQuad[1:len(ElemQuad)]
    ElemPyrd=ElemPyrd[1:len(ElemPyrd)]
    
    return ElemQuad,ElemPyrd,eleCount

def DefineElem2D90(Node,shiftX,shiftY):
    A=spatial.cKDTree(Node[:,1:3])
    
    # Find nearest
    XYPnt=np.array([0.0,0.0])
    ElemQuad=np.array([[0,0,0,0,0]])
    eleCount=1
    
    for i in range(0,len(Node)):
        idx=np.ones([1,4])*-1
        XYPnt=Node[i,1:3]
        distance1,idx1 = A.query([XYPnt[0],XYPnt[1]],k=1,distance_upper_bound=2)
        distance2,idx2 = A.query([XYPnt[0]+shiftX,XYPnt[1]],k=1,distance_upper_bound=2)
        distance3,idx3 = A.query([XYPnt[0]+shiftX,XYPnt[1]+shiftY],k=1,distance_upper_bound=2)
        distance4,idx4 = A.query([XYPnt[0],XYPnt[1]+shiftY],k=1,distance_upper_bound=2)
        idx=[idx1,idx2,idx3,idx4]
        idxTmp=np.unique(idx)

        if len(idxTmp)==4:
            ElemQuad=np.append(ElemQuad,[[eleCount,Node[idx[0],0],Node[idx[1],0],Node[idx[2],0],Node[idx[3],0] ]],axis=0)
            eleCount=eleCount+1          
    
    ElemQuad=ElemQuad[1:len(ElemQuad)]
    
    return ElemQuad,eleCount

def NodeGen3D(Node,specLenZ1):
    jmp=10**(np.ceil(np.log10(np.abs(max(Node[:,0]) + 1))))
    
    # Creating 3D Node points
    Node3D=Node
    for i in range(1,len(specLenZ1)):
        NodeTmp=np.ones(np.shape(Node))
        NodeTmp[:,0]=Node[:,0]+np.ones(np.shape(Node[:,0]))*i*jmp
        NodeTmp[:,1:3]=Node[:,1:3]
        NodeTmp[:,3]=specLenZ1[i]
        Node3D=np.append(Node3D,NodeTmp,axis=0)
        
    return Node3D,jmp

#def NodeGen3DV2(Node,zt,maxNode):
#    jmp=10**(np.ceil(np.log10(np.abs(maxNode + 1))))
#    
#    # Creating 3D Node points
#    Node3D=Node
#    NodeTmp=np.ones(np.shape(Node))
#    NodeTmp[:,0]=Node[:,0]+np.ones(np.shape(Node[:,0]))*jmp
#    NodeTmp[:,1:3]=Node[:,1:3]
#    NodeTmp[:,3]=zt
#    Node3D=np.append(Node3D,NodeTmp,axis=0)
#        
#    return Node3D,jmp

def DefineElem3D(ElemQuad,ElemPyrd,jmpNode,specLenZ1,plyStack):
    # Creating 3D pyramid elements points - 1st ply
    EleTmp=ElemPyrd[:,1:len(ElemPyrd)]
    EleTmp=EleTmp+np.ones(np.shape(ElemPyrd[:,1:len(ElemPyrd)]))*jmpNode
    ElemPyrd3D=np.append(ElemPyrd,EleTmp,axis=1)
    ElemPyrd3DPly=ElemPyrd3D
    # Generate dummy initial interface
    ElemPyrd3DInt=np.zeros(np.shape(ElemPyrd3DPly[0,:]))
    ElemPyrd3DInt=ElemPyrd3DInt.reshape(1,len(ElemPyrd3DInt))
    # Generate dummy initial interface CAM
    ElemPyrd3DCzm=np.zeros(np.shape(ElemPyrd3DPly[0,:]))
    ElemPyrd3DCzm=ElemPyrd3DCzm.reshape(1,len(ElemPyrd3DCzm))
    
    # Creating 3D quad elements points -  1st ply
    EleTmp=ElemQuad[:,1:len(ElemQuad)]
    EleTmp=EleTmp+np.ones(np.shape(ElemQuad[:,1:len(ElemQuad)]))*jmpNode
    ElemQuad3D=np.append(ElemQuad,EleTmp,axis=1)
    ElemQuad3DPly=ElemQuad3D
    # Generate dummy initial interface    
    ElemQuad3DInt=np.zeros(np.shape(ElemQuad3DPly[0,:]))   
    ElemQuad3DInt=ElemQuad3DInt.reshape(1,len(ElemQuad3DInt))
    # Generate dummy initial interface CZM 
    ElemQuad3DCzm=np.zeros(np.shape(ElemQuad3DPly[0,:]))   
    ElemQuad3DCzm=ElemQuad3DCzm.reshape(1,len(ElemQuad3DCzm))    
    
    ElemSetPly=[]
    ElemSetInt=[]
    ElemSetCzm=[]
    
    jmpElem=10**(np.ceil(np.log10(np.abs(max(ElemQuad3D[:,0]) + 1))))
    
    for i in range(1,len(specLenZ1)):
        ElemSet=[]
        # Pyramid elements
        EleTmpNds=ElemPyrd3D[:,1:len(ElemPyrd3D)]
        EleTmpNds=EleTmpNds+np.ones(np.shape(ElemPyrd3D[:,1:len(ElemPyrd3D)]))*(i-1)*jmpNode
        EleTmpNums=ElemPyrd3D[:,0]
        EleTmpNums=EleTmpNums+np.ones(np.shape(ElemPyrd3D[:,0]))*(i-1)*jmpElem
        EleTmpNums=EleTmpNums.reshape(len(EleTmpNums),1)
        EleTmpAdd=np.append(EleTmpNums,EleTmpNds,axis=1)
        if plyStack[i-1]==-1:
            ElemPyrd3DInt=np.append(ElemPyrd3DInt,EleTmpAdd,axis=0)
            ElemSetInt=np.append(ElemSetInt,ElemPyrd3DInt[:,0])
        elif plyStack[i-1]==-2:
            ElemPyrd3DCzm=np.append(ElemPyrd3DCzm,EleTmpAdd,axis=0) 
            ElemSetCzm=np.append(ElemSetCzm,ElemPyrd3DCzm[:,0])
        else:
            ElemPyrd3DPly=np.append(ElemPyrd3DPly,EleTmpAdd,axis=0)
            ElemSetPly=np.append(ElemSetPly,ElemPyrd3DPly[:,0])
        ElemSet=np.append(ElemSet,EleTmpAdd[:,0])
            
        # Quad element
        EleTmpNds=ElemQuad3D[:,1:len(ElemQuad3D)]
        EleTmpNds=EleTmpNds+np.ones(np.shape(ElemQuad3D[:,1:len(ElemQuad3D)]))*(i-1)*jmpNode
        EleTmpNums=ElemQuad3D[:,0]
        EleTmpNums=EleTmpNums+np.ones(np.shape(ElemQuad3D[:,0]))*(i-1)*jmpElem
        EleTmpNums=EleTmpNums.reshape(len(EleTmpNums),1)
        EleTmpAdd=np.append(EleTmpNums,EleTmpNds,axis=1)
        if plyStack[i-1]==-1:
            ElemQuad3DInt=np.append(ElemQuad3DInt,EleTmpAdd,axis=0)
            ElemSetInt=np.append(ElemSetInt,ElemQuad3DInt[:,0])            
        elif plyStack[i-1]==-2:
            ElemQuad3DCzm=np.append(ElemQuad3DCzm,EleTmpAdd,axis=0)
            ElemSetCzm=np.append(ElemSetCzm,ElemQuad3DCzm[:,0])                 
        else:    
            ElemQuad3DPly=np.append(ElemQuad3DPly,EleTmpAdd,axis=0)
            ElemSetPly=np.append(ElemSetPly,ElemQuad3DPly[:,0])            
        ElemSet=np.append(ElemSet,EleTmpAdd[:,0])
        writeEleSetV2(ElemSet,i)
        writeSecOriV2(ElemSet,plyStack[i-1],i)
        
    # Delete initial row
    ElemPyrd3DInt=ElemPyrd3DInt[1:len(ElemPyrd3DInt)]
    ElemQuad3DInt=ElemQuad3DInt[1:len(ElemQuad3DInt)]
    ElemPyrd3DCzm=ElemPyrd3DCzm[1:len(ElemPyrd3DCzm)]
    ElemQuad3DCzm=ElemQuad3DCzm[1:len(ElemQuad3DCzm)]   
    
    return ElemPyrd3DPly,ElemQuad3DPly,ElemPyrd3DInt,ElemQuad3DInt,ElemPyrd3DCzm,ElemQuad3DCzm,ElemSetPly,ElemSetInt,ElemSetCzm

#def DefineElem3DV2(ElemQuad,ElemPyrd,jmpNode):
#    # Creating 3D pyramid elements points - 1st ply
#    EleTmp=ElemPyrd[:,1:len(ElemPyrd)]
#    EleTmp=EleTmp+np.ones(np.shape(ElemPyrd[:,1:len(ElemPyrd)]))*jmpNode
#    ElemPyrd3D=np.append(ElemPyrd,EleTmp,axis=1)
#    ElemPyrd3D=ElemPyrd3D
#    
#    # Creating 3D quad elements points -  1st ply
#    EleTmp=ElemQuad[:,1:len(ElemQuad)]
#    EleTmp=EleTmp+np.ones(np.shape(ElemQuad[:,1:len(ElemQuad)]))*jmpNode
#    ElemQuad3D=np.append(ElemQuad,EleTmp,axis=1)
#    ElemQuad3D=ElemQuad3D
#    
#    # initialize element set
#    ElemSet=[]
#    # increment in element number
#    jmpElem=10**(np.ceil(np.log10(np.abs(max(ElemQuad3D[:,0]) + 1))))
#
#    # Pyramid elements
#    EleTmpNds=ElemPyrd3D[:,1:len(ElemPyrd3D)]
#    EleTmpNds=EleTmpNds+np.ones(np.shape(ElemPyrd3D[:,1:len(ElemPyrd3D)]))*jmpNode
#    EleTmpNums=ElemPyrd3D[:,0]
#    EleTmpNums=EleTmpNums+np.ones(np.shape(ElemPyrd3D[:,0]))*jmpElem
#    EleTmpNums=EleTmpNums.reshape(len(EleTmpNums),1)
#    EleTmpAdd=np.append(EleTmpNums,EleTmpNds,axis=1)
#    ElemPyrd3D=np.append(ElemPyrd3D,EleTmpAdd,axis=0)
#    ElemSet=np.append(ElemSet,ElemPyrd3D[:,0])
#    
#    # Quad element
#    EleTmpNds=ElemQuad3D[:,1:len(ElemQuad3D)]
#    EleTmpNds=EleTmpNds+np.ones(np.shape(ElemQuad3D[:,1:len(ElemQuad3D)]))*jmpNode
#    EleTmpNums=ElemQuad3D[:,0]
#    EleTmpNums=EleTmpNums+np.ones(np.shape(ElemQuad3D[:,0]))*jmpElem
#    EleTmpNums=EleTmpNums.reshape(len(EleTmpNums),1)
#    EleTmpAdd=np.append(EleTmpNums,EleTmpNds,axis=1)
#    ElemQuad3D=np.append(ElemQuad3D,EleTmpAdd,axis=0)
#    ElemSet=np.append(ElemSet,ElemQuad3D[:,0])            
#  
#    return ElemPyrd3D,ElemQuad3D,ElemSet

def DefineThk(specLenZ0,PlyStack,thkPly,thkInt,thkCzm):
    specLenZ1=np.array([specLenZ0])
    thk=specLenZ0
    for i in range(0,len(PlyStack)):
        if PlyStack[i]==-1:
            thk=thk+thkInt
        elif PlyStack[i]==-2:
            thk=thk+thkCzm            
        else:
            thk=thk+thkPly
        specLenZ1=np.append(specLenZ1,thk)
    
    return specLenZ1

#def writeEleSet(ElemSet):
#    f = open('EleSetFile.inp', 'w')
#    for i in range(0,len(ElemSet)):
#        elemTmp='*ELSET, GENERATE, ELSET=SET'+str(int(ElemSet[i,0]))
#        f.write("%s\n" % elemTmp)  #
#        elemTmp=str(int(ElemSet[i,1]))+','+str(int(ElemSet[i,2]))+str(int(ElemSet[i,1]))+','+str(int(ElemSet[i,2]))+str(int(ElemSet[i,1]))+','+str(int(ElemSet[i,2]))+str(int(ElemSet[i,1]))+','+str(int(ElemSet[i,2]))
#        f.write("%s\n" % elemTmp)
#    f.close() 

def writeEleSetV2(ElemSet,idt):
    ElemSet=ElemSet.astype(int)
    f = open('EleSetFile.inp', 'a+')
    elemTmp='*ELSET, GENERATE, ELSET=SET'+str(idt)
    f.write("%s\n" % elemTmp)  #
    f.close()    
    ElemSetTmp1=ElemSet[0:len(ElemSet)//8*8].reshape(len(ElemSet)//8,8)   
    with open("EleSetFile.inp", "a") as f:
        writer = csv.writer(f)     
        writer.writerows(ElemSetTmp1)
    f.close()        
    if len(ElemSet)%8>0:
        ElemSetTmp2=ElemSet[len(ElemSet)//8*8:len(ElemSet)]  
        with open("EleSetFile.inp", "a") as f:
            writer = csv.writer(f)          
            writer.writerow(ElemSetTmp2)
    f.close()
        

#def writeSecOri(ElemSet,PlyStack):
#    f = open('SecOri.inp', 'w+')
#    for i in range(0,len(ElemSet)):
#        if PlyStack[i]==-1:
#            txtTmp1='*Orientation, name=PlyOri-'+str(int(ElemSet[i,0]))
#            txtTmp2='1., 0., 0., 0., 1., 0.,'
#            txtTmp3='3, 0'
#            txtTmp4='*Solid Section, elset=SET'+str(int(ElemSet[i,0]))+', orientation=PlyOri-'+str(int(ElemSet[i,0]))+', material=matInt'
#        elif PlyStack[i]==-2:
#            txtTmp1='*Orientation, name=PlyOri-'+str(int(ElemSet[i,0]))
#            txtTmp2='1., 0., 0., 0., 1., 0.,'
#            txtTmp3='3, 0'
#            txtTmp4='*Solid Section, elset=SET'+str(int(ElemSet[i,0]))+', orientation=PlyOri-'+str(int(ElemSet[i,0]))+', material=matCzm'            
#        else:
#            txtTmp1='*Orientation, name=PlyOri-'+str(int(ElemSet[i,0]))
#            txtTmp2='1., 0., 0., 0., 1., 0.,'
#            txtTmp3='3,'+ str(PlyStack[i])
#            txtTmp4='*Solid Section, elset=SET'+str(int(ElemSet[i,0]))+', orientation=PlyOri-'+str(int(ElemSet[i,0]))+', material=matLamina'            
#        txtTmp5=','   
#         
#        f.write("%s\n" % txtTmp1)  #
#        f.write("%s\n" % txtTmp2)  #
#        f.write("%s\n" % txtTmp3)  #
#        f.write("%s\n" % txtTmp4)  #        
#        f.write("%s\n" % txtTmp5)  #        
#    f.close()   

def writeSecOriV2(ElemSet,PlyStack,idt):
    f = open('SecOri.inp', 'a+')
    if PlyStack==-1:
        txtTmp1='*Orientation, name=PlyOri-'+str(idt)
        txtTmp2='1., 0., 0., 0., 1., 0.,'
        txtTmp3='3, 0'
        txtTmp4='*Solid Section, elset=SET'+str(idt)+', orientation=PlyOri-'+str(idt)+', material=matInt'
    elif PlyStack==-2:
        txtTmp1='*Orientation, name=PlyOri-'+str(idt)
        txtTmp2='1., 0., 0., 0., 1., 0.,'
        txtTmp3='3, 0'
        txtTmp4='*Solid Section, elset=SET'+str(idt)+', orientation=PlyOri-'+str(idt)+', material=matCzm'            
    else:
        txtTmp1='*Orientation, name=PlyOri-'+str(idt)
        txtTmp2='1., 0., 0., 0., 1., 0.,'
        txtTmp3='3,'+ str(PlyStack)
        txtTmp4='*Solid Section, elset=SET'+str(idt)+', orientation=PlyOri-'+str(idt)+', material=matLamina'            
    txtTmp5=','   
     
    f.write("%s\n" % txtTmp1)  #
    f.write("%s\n" % txtTmp2)  #
    f.write("%s\n" % txtTmp3)  #
    f.write("%s\n" % txtTmp4)  #        
    f.write("%s\n" % txtTmp5)  #        
    f.close() 

def plotElem(Elem,Node):
    Elem=Elem.astype(int)
    for i in range(0,len(Elem)):
        size=len(Elem[i])
        x=[]
        y=[]
        for k in range(1,size):
            x=np.append(x,Node[Node[:,0]==Elem[i,k],1], axis=0)
            y=np.append(y,Node[Node[:,0]==Elem[i,k],2], axis=0)
        #plt.scatter(x,y)
        if size==4:
            plt.plot([x[0],x[1]],[y[0],y[1]],'r')
            plt.plot([x[1],x[2]],[y[1],y[2]],'g')
            plt.plot([x[2],x[0]],[y[2],y[0]],'k')            
        else:
            plt.plot([x[0],x[1]],[y[0],y[1]],'r')
            plt.plot([x[1],x[2]],[y[1],y[2]],'g')
            plt.plot([x[2],x[3]],[y[2],y[3]],'k')
            plt.plot([x[3],x[0]],[y[3],y[0]],'b')            
                
###############################################################################

# Inputs
#plyStack=[45,-1,-45,-1,45,-1,-45,-45,-1,45,-1,-45,-1,45]
plyStack=[45,-2,-1,-2,-45,-45,-2,-1,-2,45]
elemLen=0.075 #0.0015*np.sqrt(2)/2; # desired element size
fiberAngle=45*np.pi/180
specLenX=4
blockLen=0.5
specLenYRatio=2
thkPly=0.0075
thkInt=0.0012
thkCzm=0.0
specLenZ0=0
meshAspectRatio=1

# Non zero start point
x0=0.0
x1=blockLen
x2=x1+specLenX
x3=x2+blockLen
y0=0
yl=y0+specLenX/specLenYRatio
z0=0

# Delete prior files
os.remove('EleSetFile.inp')
os.remove('SecOri.inp')

# Derived parameters
elemLenX=np.sqrt(2)*elemLen # desired element length X
numElemX=int(np.ceil(specLenX/elemLenX)+1); #claculate desired element lentght
nodePartTemp=np.linspace(x1,x2,numElemX) # parition based on number of element requested
elemLenX=nodePartTemp[1]-nodePartTemp[0] # actual element lentgth from partitioning

elemLenY=elemLenX*meshAspectRatio; # desired element length Y
numElemY=int(np.ceil(yl/elemLenY)+1);
nodePartTemp=np.linspace(y0,yl,numElemY) # parition based on number of element requested
elemLenY=nodePartTemp[1]-nodePartTemp[0] # actual element lentgth from partitioning
specLenZ1=DefineThk(specLenZ0,plyStack,thkPly,thkInt,thkCzm)

# Shift distance
shiftX=elemLenX*np.cos(fiberAngle)**2;
shiftY=elemLenY*np.cos(fiberAngle)*np.sin(fiberAngle);

# Generate node array
NodeL, elemLenXBL,elemLenYBL=NodeGen2DV90(x0,x1,y0,yl,z0,elemLenX,elemLenY,4,numElemY)
maxNodeNum=np.max(NodeL[:,0])
NodeC=NodeGen2DV45(x1,x2,y0,yl,z0,elemLenX,elemLenY,numElemX,numElemY,shiftX,shiftY)
NodeC[:,0]=NodeC[:,0]+maxNodeNum
maxNodeNum=np.max(NodeC[:,0])
NodeR, elemLenXBR,elemLenYBR=NodeGen2DV90(x2,x3,y0,yl,z0,elemLenX,elemLenY,4,numElemY)
NodeR[:,0]=NodeR[:,0]+maxNodeNum
maxNodeNum=np.max(NodeR[:,0])

# Find boundaries
# <NCorners,NEdgeY0,NEdgeY1,NEdgeX0,NEdgeX1,NBoundary>
NCornersL,NEdgeY0L,NEdgeY1L,NEdgeX0L,NEdgeX1L,NBoundaryL=FindBoundariesV2(NodeL,x0,x1,y0,yl,numElemX,numElemY)
NCorners,NEdgeY0,NEdgeY1,NEdgeX0,NEdgeX1,NBoundary=FindBoundariesV2(NodeC,x1,x2,y0,yl,numElemX,numElemY)
NCornersR,NEdgeY0R,NEdgeY1R,NEdgeX0R,NEdgeX1R,NBoundaryR=FindBoundariesV2(NodeR,x2,x3,y0,yl,numElemX,numElemY)
#
for i in range(0,len(NEdgeX1L)):
    NodeC[(NodeC[:,0]==NEdgeX0[i,0]).nonzero()[0][0],0]=NEdgeX1L[i,0]
for i in range(0,len(NEdgeX1)):
    NodeR[(NodeR[:,0]==NEdgeX0R[i,0]).nonzero()[0][0],0]=NEdgeX1[i,0]

# Define 2D elements
# <ElemQuadL,ElemPyrdL,maxElemNoL>
ElemQuadL,maxElemNoL=DefineElem2D90(NodeL,elemLenXBL,elemLenYBL)
maxNodeNum=np.max(ElemQuadL[:,0])

ElemQuadC,ElemPyrdC,maxElemNoC=DefineElem2D45(NodeC,NEdgeY0,NEdgeY1,NEdgeX0,NEdgeX1,shiftX,shiftY)
ElemPyrdC[:,0]=ElemPyrdC[:,0]+maxNodeNum
ElemQuadC[:,0]=ElemQuadC[:,0]+maxNodeNum
maxNodeNum=np.max(ElemQuadC[:,0])

ElemQuadR,maxElemNoR=DefineElem2D90(NodeR,elemLenXBR,elemLenYBR)
ElemQuadR[:,0]=ElemQuadR[:,0]+maxNodeNum
maxNodeNum=np.max(ElemQuadR[:,0])

#
#plt.scatter(NodeL[:,1],NodeL[:,2])
plotElem(ElemQuadL,NodeL)

#plt.scatter(NodeR[:,1],NodeR[:,2])
plotElem(ElemQuadR,NodeR)

#plt.scatter(NodeC[:,1],NodeC[:,2])
plotElem(ElemQuadC,NodeC)
plotElem(ElemPyrdC,NodeC)

# Collect Nodes
Node=NodeL
Node=np.append(Node,NodeC,axis=0)
Node=np.append(Node,NodeR,axis=0)
Tmp,NodeIdx=np.unique(Node[:,0],return_index=True)
Node=Node[NodeIdx]

# Collect elements
# Quad elements
ElemQuad=ElemQuadL
ElemQuad=np.append(ElemQuad,ElemQuadC,axis=0)
ElemQuad=np.append(ElemQuad,ElemQuadR,axis=0)
# Pyramid elements
ElemPyrd=ElemPyrdC

# Find boundaries
NCorners,NEdgeY0,NEdgeY1,NEdgeX0,NEdgeX1,NBoundary=FindBoundaries(Node,x3,yl,elemLenX,elemLenY)

# Generate 3D nodes using thickness sweep
Node3D,jmpNode=NodeGen3D(Node,specLenZ1)

# Find boundaries
NCorners3D,NEdgeY03D,NEdgeY13D,NEdgeX03D,NEdgeX13D,NBoundary3D=FindBoundaries(Node3D,x3,yl,elemLenX,elemLenY)

# Define 3D elements
ElemPyrd3DPly,ElemQuad3DPly,ElemPyrd3DInt,ElemQuad3DInt,ElemPyrd3DCzm,ElemQuad3DCzm,ElemSetPly,ElemSetInt,ElemSetCzm=DefineElem3D(ElemQuad,ElemPyrd,jmpNode,specLenZ1,plyStack)

## Write data to csv file.
np.savetxt("Node3D.csv", Node3D, delimiter=",", fmt=('%1.2i','%1.6f','%1.6f','%1.6f'))
#
np.savetxt("ElemPyrd3DPly.csv", ElemPyrd3DPly, delimiter=",", fmt='%1.2i')
np.savetxt("ElemQuad3DPly.csv", ElemQuad3DPly, delimiter=",", fmt='%1.2i')
#if min(plyStack)==-1:
np.savetxt("ElemPyrd3DInt.csv", ElemPyrd3DInt, delimiter=",", fmt='%1.2i')
np.savetxt("ElemQuad3DInt.csv", ElemQuad3DInt, delimiter=",", fmt='%1.2i')
#if min(plyStack)==-1:
np.savetxt("ElemPyrd3DCzm.csv", ElemPyrd3DCzm, delimiter=",", fmt='%1.2i')
np.savetxt("ElemQuad3DCzm.csv", ElemQuad3DCzm, delimiter=",", fmt='%1.2i')
#
## Print stats
#print('Total nodes: ',max(Node3D[:,0]))
#print('Total pyramid elements - Ply: ',max(ElemPyrd3DPly[:,0]))
#print('Total quad elements - Ply: ',max(ElemQuad3DPly[:,0]))
##if min(plyStack)==-1:
#print('Total pyramid elements - Interface: ',max(ElemPyrd3DInt[:,0]))
#print('Total quad elements - Interface: ',max(ElemQuad3DInt[:,0]))
##if min(plyStack)==-1:
#print('Total pyramid elements - Czm: ',max(ElemPyrd3DCzm[:,0]))
#print('Total quad elements - Czm: ',max(ElemQuad3DCzm[:,0]))
#print('Total elements: ',max(ElemPyrd3DPly[:,0])+max(ElemQuad3DPly[:,0]))
