'''
Created on Jun 29, 2013

@author: Wenhao Sun
'''

import fractions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
import operator
import numpy as np
import numpy.linalg as npl
from numpy import pi


"""
Dear Surface scientist, 

Thank you for taking interest in my algorithm. Below is a working implementation of the 3-atom screening heuristic. 

There are many potential modifications and extensions which could make this algorithm more powerful.
Here, it is provided in a basic form, but it captures the physical essence of the 3-atom heuristic.    

This was written prior to the implementation of the pymatgen.core.surface package. 
There may be implementations of existing python methods that can simplify this code. 

Best regards,
Wenhao Sun
October, 2017 
"""


class ScreenSurf():
    def __init__(self,bulkstruct):
        
        RefinedStructure=SpacegroupAnalyzer(bulkstruct,0.03).get_refined_structure()
        if StructureMatcher().fit(bulkstruct,RefinedStructure):
            bulkstruct=RefinedStructure
            self.bulkstruct=self.getconventional(bulkstruct)
        else:
            self.bulkstruct=bulkstruct
        
    def ScreenSurfaces(self,natomsinsphere=50,keep=1,samespeconly=True,ignore=[]):
        '''
        Screen for low energy surfaces with 3-atom histogram heuristic
        
        ScreenSurfaces uses the following parameters:
            natomsinsphere = Increases radius of sphere until there are #natomsinsphere in the sphere. For simple structures this can be as low as 20. 
            keep = What percentage of the histogram to keep up to. Using keep=0.8 means we only keep the most common 80% of planes. 
            samespeconly = This is whether or not to search for only atoms of the same element in the sphere. 
            ignore = these are atoms to ignore. These often include small atoms on complex anions (e.g., O on CO_3^2-), or other moieties, like hydrogen. 
        '''
        print "\n*Running Screening Heuristic"
        liststr,longlist=self.planesearch(self.bulkstruct,natomsinsphere=natomsinsphere,keep=keep,samespeconly=samespeconly,ignore=ignore)
        print liststr
        if len(liststr)==0:
            print 'No Integer Miller Indices'
        else:
            index=np.array(self.indexfromstr(liststr[0][0]))
        liststr.pop(0)
        for str in liststr:
            ind=np.array(self.indexfromstr(str[0]))
            index=np.vstack((index,ind))
        return index
    
    def planesearch(self,bulk,natomsinsphere=50,keep=0.8,samespeconly=False,ignore=[],maxmult=9):
        '''
        People typically assume the 'low-index' surfaces are low in energy. However, since the 
        definition of the unit cell is a *human* (arguably arbitrary) convention, not a natural convention, there ought to 
        exist a unit-cell-independent approach to identify low-energy surface cuts. 
        
        The hypothesis here is that low-energy surfaces tend to be planes such that there is a high concentration
        of atoms that lie upon plane. The algorithm here is to enumerate all combinations of three atoms, enumerating from the 
        *atoms* with the highest symmetry within the unit cell. There is a unit-cell-independent way to do this

        1. We identify all non-equivalent atoms in the unit cell. 
        2. We pick the atoms with the highest symmetry
        3. We find planes with all other atoms of the same chemical identity, choosing atoms within a radius that allows proper sampling. 
        4. We build a histogram of index occurrences, assuming that indices with high occurrence have higher density of atoms,
           and are thus lower in energy to break. 
           
        The domination of low-energy planes is thus that which coincides with the planes with the highest density of atoms, 
        and therefore the planes with the least number of bonds broken.   
        '''
        
        basis=bulk.lattice.matrix
        convstruct=self.getconventional(bulk)
        '1. We identify all non-equivalent atoms in the unit cell.'
        SS=SpacegroupAnalyzer(bulk).get_symmetrized_structure()
        uniqsites=[]
        for eq in SS.equivalent_sites:
            if str(eq[0].specie) not in ignore:
                uniqsites.append(eq[0])
        
        '2. Find an *r* that gives some minimum number of atoms in the sphere for each atom'
        'Start with r=4, increase until natomsinsphere is achieved'
        r=4
        for site in uniqsites:
            pt=site.coords
            insphere=bulk.get_sites_in_sphere(pt, r, True)
            while len(insphere)<=natomsinsphere: 
                r=r+0.1
                insphere=bulk.get_sites_in_sphere(pt, r, True)
        
        print "Radius of sphere for screening heuristic = "+str(r)
        print "Fraction of planes kept:" +str(keep)
        print "Search on same species only?: "+str(samespeconly)
        print "Ignoring elements?: "+str(ignore)
        
        '''3. We find planes with all other atoms, either of the same chemical identity or not, depending on
           samespeconly variable. We choose atoms within a radius that is some 
           large amount greater than the size of the unit cell'''
        
        totaldict=dict()
        for site in uniqsites:
            millercheck=dict()
            symfam=dict()
            sitecollection=[]
            pt=site.coords
            insphere=bulk.get_sites_in_sphere(pt, r, True)
            for x in insphere:
                if samespeconly:
                    if x[0].specie==site.specie:
                        sitecollection.append(x[0].coords)
                else: 
                    if str(x[0].specie) not in ignore:
                        sitecollection.append(x[0].coords)
            
            #Generate all combinations of three atoms in this list.         
            M=len(sitecollection)
            POS=[]
            for ii in range(1,M+1):
                for jj in range (ii+1,M+1):
                    for kk in range(jj+1,M+1):
                        ncr=[ii ,jj ,kk]
                        POS.append(np.array(ncr))
                        
            for ii in range(0,len(POS)):
                Setof3Atoms=np.array([sitecollection[POS[ii][0]-1], sitecollection[POS[ii][1]-1], sitecollection[POS[ii][2]-1]])
                
                for aa in range(0,3):
                    v1=Setof3Atoms[np.mod(aa+1,3)]-Setof3Atoms[np.mod(aa,3)]
                    v2=Setof3Atoms[np.mod(aa+2,3)]-Setof3Atoms[np.mod(aa,3)]
                    T=self.ang(v1, v2)
                    if T <= pi/2:
                        bb=aa
                        break
                v1xv2=np.cross(v1,v2)
                
                if npl.norm(v1xv2)==0:
                    continue
                if self.ang(v1,v2)<1E-3:
                    continue
                if v1xv2[0] < 0:
                    vtemp=v1
                    v1=v2
                    v2=vtemp
                    v1xv2=v1xv2*-1
                elif v1xv2[0]==0:
                    if v1xv2[1]<0:
                        vtemp=v1
                        v1=v2
                        v2=vtemp
                        v1xv2=v1xv2*-1
                    elif v1xv2[1]==0:
                        if v1xv2[2]<0:
                            vtemp=v1
                            v1=v2
                            v2=vtemp
                            v1xv2=v1xv2*-1
                else:
                    continue
                
                millerindex=self.getmillerfrom2v(basis, v1, v2,False)
                if millerindex==None or (millerindex > 10).any():
                    continue
                found=False
                strm=self.indextostr(millerindex)
                for fam in symfam.keys():
                    if strm in symfam[fam]:
                        strm=fam
                        millerindex=self.indexfromstr(strm)
                        found=True
                if found==False:              
                    C=self.spggetsymfam(convstruct,millerindex)
                    millerindex=C[-1]
                    strm=self.indextostr(millerindex)
                    strsymfam=[]
                    for x in C:
                        strsymfam.append(self.indextostr(x))
                    symfam[strm]=strsymfam
                    
                duplicate=0
                replace=0
                strI=self.indextostr(millerindex)
                
                if strI in millercheck.keys():
                    millercheck[strI]=millercheck[strI]+1
                else:
                    millercheck[strI]=1
            
            sorted_millercheck = sorted(millercheck.iteritems(), key=operator.itemgetter(1),reverse=True)
            for miller in millercheck.keys():
                if miller in totaldict.keys():
                    totaldict[miller]=totaldict[miller]+millercheck[miller]
                else:
                    totaldict[miller]=millercheck[miller]
        sorttotal=sorted(totaldict.iteritems(), key=operator.itemgetter(1),reverse=True)
        
        total=0
        for ind in sorttotal:
            total=total+ind[1]
        keeplist=[]
        cdf=0
        alltotal=0
        for ind in sorttotal:
            alltotal=alltotal+ind[1]+0.0
        
        for ind in sorttotal:
            cdf=cdf+ind[1]
            #In this section, we limit the length of the indices. 
            #Don't allow indices greater with more than 1 digit 
            short=True
            mult=1
            for ii in self.indexfromstr(ind[0]):
                mult=mult*ii
            for aa in self.indexfromstr(ind[0]):
                if len(str(aa).replace("-",""))>1:
                    short=False
            if cdf < total*keep:
                if short==True and mult<=maxmult:
                    keeplist.append((ind[0],'{percent:.1%}'.format(percent=ind[1]/alltotal)))
            else:
                break
        return keeplist,sorttotal
    
    def indexfromstr(self,str):
        '''Follows format of str: [h,k,l]'''
        index=str.replace("("," ").replace(")"," ").replace(","," ")
        if self._representsint(index.split()[0]) and self._representsint(index.split()[1]) and  self._representsint(index.split()[2]): 
            h=int(index.split()[0])
            k=int(index.split()[1])
            l=int(index.split()[2])
            return np.array([h,k,l]) 
        else:
            print str
            raise SystemError("str is not in index form - not integers")
    
    def getconventional(self,bulkstruct):
        return SpacegroupAnalyzer(bulkstruct,0.03).get_conventional_standard_structure()
    
    
    def ang(self, v1, v2,acute=False):
        if npl.norm(v1)==0 or npl.norm(v2)==0:
            raise StandardError('One of your vectors has length zero')
        x=np.dot(v1,v2)/(npl.norm(v1)*npl.norm(v2))
        a= np.float('%.12f' %x)
        angle=np.arccos(a)
        if acute==True:
            if angle>pi/2:
                angle=pi-angle
        return angle
    
    def getmillerfrom2v(self,basis,v1,v2,errfrac=True):
        
        TM=np.eye(3,3)
        millerv1=np.dot(TM,np.transpose(v1))
        millerv2=np.dot(TM,np.transpose(v2))
        
        if self.ang(millerv1,millerv2,True)<1E-2:
            return None
        
        millerv1=np.dot(npl.inv(np.transpose(basis)),millerv1)
        millerv2=np.dot(npl.inv(np.transpose(basis)),millerv2)
        
        millerindex=np.transpose(np.cross(millerv1,millerv2))
        if np.isnan(millerindex).any():
                print millerindex
                print v1,millerv1
                print v2,millerv2
                raise SystemError
        
        md=np.zeros((3,1))
        for ni in range(0,3):
            md[ni]=fractions.Fraction(millerindex[ni]).limit_denominator(12).denominator
        
        MLCM=1
        for ii in md:
            MLCM=MLCM*ii / fractions.gcd(MLCM,ii)
        
        millerindex=millerindex*MLCM;
        roundedmillerindex=np.array([ round(elem,1) for elem in millerindex])
        diffmill=roundedmillerindex-millerindex
        if npl.norm(diffmill)>0.005:
            if errfrac==True:
                raise SystemError('Error: Miller Index is Fractional')
            else:
                return None
        
        for elem in range(0,3):
            if roundedmillerindex[elem]==-0:
                roundedmillerindex[elem]=0
        
        millerindex=roundedmillerindex
        MGCD=np.abs(fractions.gcd(fractions.gcd(millerindex[0],millerindex[1]),millerindex[2]))
        
        if MGCD==0:
            print "Zero in the Miller greatest common denominator"
            print millerindex
            print millerv1, millerv2
            return None
        
        millerindex=millerindex/MGCD
        return millerindex    
    
    def indextostr(self,I):
        try:
            return "("+str(int(I[0]))+","+str(int(I[1]))+","+str(int(I[2]))+")"
        except:
            print I
            raise SystemError("Fail Index")
    
    def spggetsymfam(self,structure, mindex):
        basis=np.array(structure.lattice.matrix)
        try:
            A,useless=self.get2vectsinplane(basis,mindex)
            v1=A[0]
            v2=A[1]
            if np.isnan(v2).any():
                print v2
                print mindex
                raise SystemError
        except:
            print "Failed get2vectsinplane"
            return mindex
        millerindex=self.getmillerfrom2v(basis, v1, v2)
        
        if (millerindex==mindex*-1).all():
            v2=A[0]
            v1=A[1]
        C=[]
        latt = structure.lattice
        symmops = SpacegroupAnalyzer(structure, 0.01).get_symmetry_operations(True)
        for op in symmops:
            for site in structure:
                newv1 = op.apply_rotation_only(v1)
                newv2 = op.apply_rotation_only(v2)
                millerindex=np.array(self.getmillerfrom2v(basis, newv1, newv2,errfrac=False))
                if millerindex==None:
                    continue
                exists=False
                for Ccheck in C:
                    c=0
                    for ind in range(0,3):
                        if Ccheck[ind]==millerindex[ind]:
                            c+=1
                    if c==3:
                        exists=True                    
                if exists==False:
                    C.append(millerindex)    
        
        retC=np.zeros((len(C),3))
        for Cx in range(0,len(C)):
            retC[Cx]=C[Cx]
        a=retC
        ind=np.lexsort((a[:,2],a[:,1],a[:,0]))
        return a[ind]
    
    def get2vectsinplane(self,basis, maxindex):
        #% The two returned vectors are in cartesian coordinates
        if len(maxindex) != 3:
            raise StandardError('Error need 3 elements in your maxindex')
    
        zeroind=np.where(maxindex==0)[0]
        numzero = len(zeroind)
        
        h = float(maxindex[0])
        k = float(maxindex[1])
        l = float(maxindex[2])
       
        if numzero == 0:
            hd=fractions.Fraction(1/float(h)).limit_denominator(12).denominator
            kd=fractions.Fraction(1/float(k)).limit_denominator(12).denominator
            ld=fractions.Fraction(1/float(l)).limit_denominator(12).denominator
            lst=[hd,kd,ld]
            multfact=1
            for ii in lst:
                multfact=multfact*ii / fractions.gcd(multfact,ii)

            p1=np.array([multfact/h,0,0])
            p2=np.array([0,multfact/k,0])
            p3=np.array([0,0,multfact/l])
            
            P=np.array([p1,p2,p3])
            for aa in range(0,3):
                v1=P[np.mod(aa+1,3),:]-P[np.mod(aa,3),:]
                v2=P[np.mod(aa+2,3),:]-P[np.mod(aa,3),:]
                T=self.ang(v1,v2);
                if T<=pi/2:
                    twovects=np.array([v1,v2]);
                    break;
            
        elif numzero==1:
            ind=list()
            P=list()
            for jj in range(0,3):
                if jj==zeroind:
                    p1=np.array([0, 0, 0])
                    p1[jj]=1
                    v1=p1
                else:
                    ind.append(jj)
                    P.append(maxindex[jj])
            
            ad=fractions.Fraction(1/float(P[0])).limit_denominator(12).denominator
            bd=fractions.Fraction(1/float(P[1])).limit_denominator(12).denominator
            
            
            lst=[ad,bd]
            multfact=1
            for ii in lst:
                multfact=multfact*ii / fractions.gcd(multfact,ii)
                
            points=np.zeros((2,3));
            for mm in range(0,2):
                pointtemp=np.array([0, 0, 0])
                pointtemp[ind[mm]]=multfact/P[mm]
                points[mm,:]=pointtemp
            v2=points[1]-points[0]
            
            
            twovects=np.array([v1,v2])
            P=np.array([p1, points[0],points[1]])
            
        elif numzero==2:
            maxindex=maxindex/npl.norm(maxindex);
            b = []
            for i in maxindex:
                b.append(abs(float(i)))
            maxindex=b        
            if maxindex == [1, 0, 0]:
                twovects=np.array([[0, 0, 1],[0, 1, 0]])
            elif maxindex==[0, 1, 0]:
                twovects=np.array([[1, 0, 0],[0, 0, 1]])
            elif maxindex==[0, 0, 1]:
                twovects=np.array([[1, 0, 0],[0, 1, 0]])

            P=np.array([twovects[0], twovects[1], [0, 0,0]]);
            
        twovects[0]=twovects[0]/fractions.gcd(fractions.gcd(twovects[0,0], twovects[0,1]),twovects[0,2])
        twovects[1]=twovects[1]/fractions.gcd(fractions.gcd(twovects[1,0], twovects[1,1]),twovects[1,2])
        twovects=np.dot(twovects,basis);
        for v in twovects:
            for x in v:
                if x > 10:
                    pass
                    #raise SystemError("Too big vectors??")
                    
        return twovects,P
    
    def _representsint(self,s):
        try: 
            int(s)
            return True
        except ValueError:
            return False
