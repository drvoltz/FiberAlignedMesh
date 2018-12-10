
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial

def NodeGen2D(specLenX,specLenZ0,elemLenX,elemLenY,shiftX,shiftY):
    # Nodal coordinated
    nodeX1=np.arange(0,specLenX+elemLenX,elemLenX);
    nodeY1=np.zeros(np.shape(nodeX1)) #np.arange(0,specLenY+elemLenY,elemLenY);
    #
    nodeX2=np.arange(shiftX,specLenX+shiftX,elemLenX);
    nodeY2=np.zeros(np.shape(nodeX2))+shiftY
    #
    # Create all nodes
    count=1;
    #jmp=10**(np.floor(np.log10(np.abs(specLenX) + 1)));
    Node=np.array([[0,0,0,0]])
    for j in range(0,int(numElemY)):
        for i in range(0,len(nodeX1)):
            Node=np.append(Node,[[count+i,nodeX1[i],nodeY1[i]+j*elemLenY,specLenZ0]],axis=0)
        count=len(Node)
        
        for i in range(0,len(nodeX2)):
            Node=np.append(Node,[[count+i,nodeX2[i],nodeY2[i]+j*elemLenY,specLenZ0]],axis=0)
        count=len(Node)      
        
    # last line
    for i in range(0,len(nodeX1)):
        Node=np.append(Node,[[count+i,nodeX1[i],nodeY1[i]+(j+1)*elemLenY,specLenZ0]],axis=0)
            
    Node=Node[1:len(Node)]
    return Node

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

def DefineElem2D(Node,NBtmEdge,NTopEdge,NLeftEdge,NRightEdge):
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
    
    ElemSet=[[ 1,min(ElemPyrd3D[:,0]),max(ElemQuad3D[:,0]) ]]
    
    jmpElem=10**(np.ceil(np.log10(np.abs(max(ElemQuad3D[:,0]) + 1))))
    
    for i in range(1,len(specLenZ1)):
        # Pyramid elements
        EleTmpNds=ElemPyrd3D[:,1:len(ElemPyrd3D)]
        EleTmpNds=EleTmpNds+np.ones(np.shape(ElemPyrd3D[:,1:len(ElemPyrd3D)]))*(i-1)*jmpNode
        EleTmpNums=ElemPyrd3D[:,0]
        EleTmpNums=EleTmpNums+np.ones(np.shape(ElemPyrd3D[:,0]))*(i-1)*jmpElem
        EleTmpNums=EleTmpNums.reshape(len(EleTmpNums),1)
        EleTmpAdd=np.append(EleTmpNums,EleTmpNds,axis=1)
        if plyStack[i-1]==-1:
            ElemPyrd3DInt=np.append(ElemPyrd3DInt,EleTmpAdd,axis=0)
        elif plyStack[i-1]==-2:
            ElemPyrd3DCzm=np.append(ElemPyrd3DCzm,EleTmpAdd,axis=0)            
        else:
            ElemPyrd3DPly=np.append(ElemPyrd3DPly,EleTmpAdd,axis=0)
            
        ElemMin=min(EleTmpAdd[:,0])
        # Quad element
        EleTmpNds=ElemQuad3D[:,1:len(ElemQuad3D)]
        EleTmpNds=EleTmpNds+np.ones(np.shape(ElemQuad3D[:,1:len(ElemQuad3D)]))*(i-1)*jmpNode
        EleTmpNums=ElemQuad3D[:,0]
        EleTmpNums=EleTmpNums+np.ones(np.shape(ElemQuad3D[:,0]))*(i-1)*jmpElem
        EleTmpNums=EleTmpNums.reshape(len(EleTmpNums),1)
        EleTmpAdd=np.append(EleTmpNums,EleTmpNds,axis=1)
        if plyStack[i-1]==-1:
            ElemQuad3DInt=np.append(ElemQuad3DInt,EleTmpAdd,axis=0)
        elif plyStack[i-1]==-2:
            ElemQuad3DCzm=np.append(ElemQuad3DCzm,EleTmpAdd,axis=0)             
        else:    
            ElemQuad3DPly=np.append(ElemQuad3DPly,EleTmpAdd,axis=0)
        ElemMax=max(EleTmpAdd[:,0])
        #
        ElemSet=np.append(ElemSet,[[i,ElemMin,ElemMax]],axis=0)

    # Delete initial row        
    ElemSet=ElemSet[1:len(ElemSet)]
    ElemPyrd3DInt=ElemPyrd3DInt[1:len(ElemPyrd3DInt)]
    ElemQuad3DInt=ElemQuad3DInt[1:len(ElemQuad3DInt)]
    ElemPyrd3DCzm=ElemPyrd3DCzm[1:len(ElemPyrd3DCzm)]
    ElemQuad3DCzm=ElemQuad3DCzm[1:len(ElemQuad3DCzm)]    
    
    return ElemPyrd3DPly,ElemQuad3DPly,ElemPyrd3DInt,ElemQuad3DInt,ElemPyrd3DCzm,ElemPyrd3DCzm,ElemSet

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

def writeEleSet(ElemSet):
    f = open('EleSetFile.inp', 'w')
    for i in range(0,len(ElemSet)):
        elemTmp='*ELSET, GENERATE, ELSET=SET'+str(int(ElemSet[i,0]))
        f.write("%s\n" % elemTmp)  #
        elemTmp=str(int(ElemSet[i,1]))+','+str(int(ElemSet[i,2]))+',1'
        f.write("%s\n" % elemTmp)
    f.close() 

def writeSecOri(ElemSet,PlyStack):
    f = open('SecOri.inp', 'w')
    for i in range(0,len(ElemSet)):
        if PlyStack[i]==-1:
            txtTmp1='*Orientation, name=PlyOri-'+str(int(ElemSet[i,0]))
            txtTmp2='1., 0., 0., 0., 1., 0.,'
            txtTmp3='3, 0'
            txtTmp4='*Solid Section, elset=SET'+str(int(ElemSet[i,0]))+', orientation=PlyOri-'+str(int(ElemSet[i,0]))+', material=matInt'
        elif PlyStack[i]==-2:
            txtTmp1='*Orientation, name=PlyOri-'+str(int(ElemSet[i,0]))
            txtTmp2='1., 0., 0., 0., 1., 0.,'
            txtTmp3='3, 0'
            txtTmp4='*Solid Section, elset=SET'+str(int(ElemSet[i,0]))+', orientation=PlyOri-'+str(int(ElemSet[i,0]))+', material=matCzm'            
        else:
            txtTmp1='*Orientation, name=PlyOri-'+str(int(ElemSet[i,0]))
            txtTmp2='1., 0., 0., 0., 1., 0.,'
            txtTmp3='3,'+ str(PlyStack[i])
            txtTmp4='*Solid Section, elset=SET'+str(int(ElemSet[i,0]))+', orientation=PlyOri-'+str(int(ElemSet[i,0]))+', material=matLamina'            
        txtTmp5=','   
         
        f.write("%s\n" % txtTmp1)  #
        f.write("%s\n" % txtTmp2)  #
        f.write("%s\n" % txtTmp3)  #
        f.write("%s\n" % txtTmp4)  #        
        f.write("%s\n" % txtTmp5)  #        
    f.close()   
    
###############################################################################

# Inputs
#plyStack=[45,-1,-45,-1,45,-1,-45,-45,-1,45,-1,-45,-1,45]
plyStack=[45,-2,-1,-2,-45,-2,-1,-2,45]
elemLen=0.0075 #0.0015*np.sqrt(2)/2; # desired element size
fiberAngle=45*np.pi/180;
specLenX=0.5;
specLenY=0.25;
thkPly=0.0075
thkInt=0.0012
thkCzm=0.0
specLenZ0=0;
#specLenZ1=[0,thkPly]
meshAspectRatio=1

# Derived parameters
elemLenX=np.sqrt(2)*elemLen
elemLenY=elemLenX*meshAspectRatio;
numElemX=np.ceil(specLenX/elemLenX);
numElemY=np.ceil(specLenY/elemLenY);
specLenZ1=DefineThk(specLenZ0,plyStack,thkPly,thkInt,thkCzm)

# Recaculate parameters
elemLenX=specLenX/numElemX;
elemLenY=specLenY/numElemY;

# Shift distance
shiftX=elemLenX*np.cos(fiberAngle)**2;
shiftY=elemLenY*np.cos(fiberAngle)*np.sin(fiberAngle);

# Generate node array
Node=NodeGen2D(specLenX,specLenZ0,elemLenX,elemLenY,shiftX,shiftY)

# Find boundaries
NCorners,NEdgeY0,NEdgeY1,NEdgeX0,NEdgeX1,NBoundary=FindBoundaries(Node,specLenX,specLenY,elemLenX,elemLenY)

# Generate 3D nodes using thickness sweep
Node3D,jmpNode=NodeGen3D(Node,specLenZ1)

# Define 2D elements
ElemQuad,ElemPyrd,maxElemNo=DefineElem2D(Node,NEdgeY0,NEdgeY1,NEdgeX0,NEdgeX1)

# Define 3D elements
ElemPyrd3DPly,ElemQuad3DPly,ElemPyrd3DInt,ElemQuad3DInt,ElemPyrd3DCzm,ElemQuad3DCzm,ElemSet=DefineElem3D(ElemQuad,ElemPyrd,jmpNode,specLenZ1,plyStack)

# Write element set
writeEleSet(ElemSet)

# Write element set
writeSecOri(ElemSet,plyStack)

# Write data to csv file.
np.savetxt("Node3D.csv", Node3D, delimiter=",", fmt=('%1.2i','%1.6f','%1.6f','%1.6f'))    
np.savetxt("ElemPyrd3DPly.csv", ElemPyrd3DPly, delimiter=",", fmt='%1.2i')
np.savetxt("ElemQuad3DPly.csv", ElemQuad3DPly, delimiter=",", fmt='%1.2i')
#if min(plyStack)==-1:
np.savetxt("ElemPyrd3DInt.csv", ElemPyrd3DInt, delimiter=",", fmt='%1.2i')
np.savetxt("ElemQuad3DInt.csv", ElemQuad3DInt, delimiter=",", fmt='%1.2i')
#if min(plyStack)==-1:
np.savetxt("ElemPyrd3DCzm.csv", ElemPyrd3DCzm, delimiter=",", fmt='%1.2i')
np.savetxt("ElemQuad3DCzm.csv", ElemQuad3DCzm, delimiter=",", fmt='%1.2i')

# Print stats
print('Total nodes: ',max(Node3D[:,0]))
print('Total pyramid elements - Ply: ',max(ElemPyrd3DPly[:,0]))
print('Total quad elements - Ply: ',max(ElemQuad3DPly[:,0]))
#if min(plyStack)==-1:
print('Total pyramid elements - Interface: ',max(ElemPyrd3DInt[:,0]))
print('Total quad elements - Interface: ',max(ElemQuad3DInt[:,0]))
#if min(plyStack)==-1:
print('Total pyramid elements - Czm: ',max(ElemPyrd3DCzm[:,0]))
print('Total quad elements - Czm: ',max(ElemQuad3DCzm[:,0]))
print('Total elements: ',max(ElemPyrd3DPly[:,0])+max(ElemQuad3DPly[:,0]))
