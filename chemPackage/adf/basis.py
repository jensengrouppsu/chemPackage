from __future__ import print_function, division
from numpy import array

def collect_basis(self, f, indices):
    '''Collects the different types of basis functions.'''

    self.basis = {}

    if 'SFO BASIS' in indices:

        # Initialize
        self.basis.update({'SFO': []})

        # Find the start and end of the SFO basis block
        s = indices['SFO BASIS']
        try:
            e = [i for i in xrange(s,min(s+5000,len(f)-1)) if f[i]==''][0]
        except IndexError:
            return
        # Get SFO basis
        try:
            aoid = [f[i].split()[9]+' '+f[i].split()[5] for i in xrange(s,e,2)]
            orbital = [f[i].split()[7]+' '+f[i].split()[8] for i in xrange(s,e,2)]
        except IndexError:
            return
        # Add values to key
        [self.basis['SFO'].append(tuple((aoid[i], orbital[i]))) for i in xrange(len(aoid))]

    if 'STO BASIS' in indices:
        # Initialize
        self.basis.update({'STO': {}})

        # Find start and end of STO basis block
        s = indices['STO BASIS']
        
        try:
            #Xing modify
            #e = [i-6 for i in xrange(s,min(s+5000,len(f)-1)) if '****' in f[i]][0] 
            e = [i for i in xrange(s,min(s+5000,len(f)-1)) if 'Total number of charge' in f[i]][0]
        except IndexError:
            return
        
        # Find the index of the start of the types of atoms block
        atindx = [i for i, x in enumerate(f[s:e], s) if f[i][0:4].strip(' ')!='']

        # Collect atom type and atom numbers 
        attype = [f[i][0:4].strip(' ') for i in atindx]
        #atnumb = [f[i].split()[1:] for i in atindx]
        atnumb=[]
        atindxs=[]
        for i in atindx:
            temp = f[i].split()[1:]
            if f[i+1][:35].strip(' ')=='-':
                atnumb.append(temp)
                atindxs.append(i)
                temp=[]
            while f[i+1][:35].strip(' ')!='-':
                i = i + 1  
                temp = temp + f[i].split()[:]    
            if len(temp) != 0: 
                atnumb.append(temp)
                atindxs.append(i)    
        atindxe = atindx
        atindxe.append(e)
        # Add a key for each atom
        [[self.basis['STO'].update({atnumb[i][j]+' '+attype[i]: {}}) for j in xrange(len(atnumb[i]))] for i in xrange(len(attype))]

        # Collect X, Y, Z, and alpha for each type
        try:
            X = [[int(f[i][10:33].split()[0]) for i in xrange(atindxs[a]+2,atindxe[a+1]) if len(f[i][10:33].split())>4] for a in range(len(attype))]
            Y = [[int(f[i][10:33].split()[1]) for i in xrange(atindxs[a]+2,atindxe[a+1]) if len(f[i][10:33].split())>4] for a in range(len(attype))]
            Z = [[int(f[i][10:33].split()[2]) for i in xrange(atindxs[a]+2,atindxe[a+1]) if len(f[i][10:33].split())>4] for a in range(len(attype))]
            R = [[int(f[i][10:33].split()[3]) for i in xrange(atindxs[a]+2,atindxe[a+1]) if len(f[i][10:33].split())>4] for a in range(len(attype))]
            alf = [[float(f[i][10:33].split()[4]) for i in xrange(atindxs[a]+2,atindxe[a+1]) if len(f[i][10:33].split())>4] for a in range(len(attype))]
        except ValueError or IndexError:
            return
        # Figure out the orbital type from the X, Y, Z, and R values
        sL = ['S', 'P', 'D', 'F']
        sX = ['', 'x', 'x2', 'x3', 'x4', 'x5']
        sY = ['', 'y', 'y2', 'y3', 'y4', 'y5']
        sZ = ['', 'z', 'z2', 'z3', 'z4', 'z5']
        obtyp = [['S' if X[i][j]+Y[i][j]+Z[i][j]==0 else sL[X[i][j]+Y[i][j]+Z[i][j]]+':'+sX[X[i][j]]+sY[Y[i][j]]+sZ[Z[i][j]]
                  for j in xrange(len(X[i]))] for i in xrange(len(X))]
        val = [[obtyp[i][:j].count(obtyp[i][j])+1 for j in xrange(len(obtyp[i]))] for i in xrange(len(obtyp))]
        obtyp = [[str(val[i][j])+' '+obtyp[i][j] for j in xrange(len(obtyp[i]))] for i in xrange(len(obtyp))]
        # Add information to each atom basis
        [[[self.basis['STO'][atnumb[i][j]+' '+attype[i]].update({obtyp[i][k]: tuple((X[i][k], Y[i][k], Z[i][k], R[i][k], alf[i][k]))})
           for k in xrange(len(X[i]))] for j in xrange(len(atnumb[i]))] for i in xrange(len(attype))]
