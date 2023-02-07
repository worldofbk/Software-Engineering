from sympy import *
#---------------------
#Copy Matrix Function
def zeromat(r,c):
  m=[]
  for i in range(c):
    temp=[0.0]*r
    m.append(temp)
  return m

def copymat(m):
    m1=[]
    r = len(m)
    c = len(m[0])
    m1=zeromat(r,c)
    for i in range(r):
        for j in range(c):
            m1[i][j] = m[i][j]
    return m1

#---------------------
#Print Matrix Function
def print_mat(mat):
  matli=[]
  for i in range(len(mat)):
    temp=[]
    for j in range(len(mat[i])):
      temp.append(len(str(mat[i][j])))
      matli.append(temp)
  maxlen=max(matli)

  print("Wronskian Matrix")
  print()
  for i in range(len(mat)):
    for j in range(len(mat[i])):
      print("|",end="")
      print(mat[i][j],(maxlen[j]-(len(str(mat[i][j]))))*" ",end="")
    print("|")  
    print()

import sympy as sp
x = sp.Symbol('x')

#-------------------------
#Wronskian Matrix Function

def wronskianmat(f): #[[x**2,x**3][2*x,3*x**2]]
  import sympy as sp
  x = sp.Symbol('x')
  num = len(f)
  mat=[]
  mat.append(f)
  for i in range(1,num):
    temp=[]
    for j in mat[i-1]:
      temp.append(sp.diff(j))
    mat.append(temp)
  print()
  return mat
#--------------
#Function Input
f=[]
print("Enter number of equations: ")
for i in range(int(input())):
  s = input("Enter an Equation : ")
  f.append(sp.sympify(s)) #converts the input into mathematical function
print("Input Equations :")

for i in range(len(f)):
  print(f[i])
  print()


#-----------------------

finalmat=(wronskianmat(f))

def det2(matdetfun):
  if(len(matdetfun)==2):
    return((sp.simplify(matdetfun[0][0])*sp.simplify(matdetfun[1][1]))-(sp.simplify(matdetfun[1][0])*sp.simplify(matdetfun[0][1])))
def det3(matdetfun):
    return ((sp.simplify(matdetfun[0][0])*((sp.simplify(matdetfun[1][1])*sp.simplify(matdetfun[2][2]))-(sp.simplify(matdetfun[2][1])*sp.simplify(matdetfun[1][2]))))-
            (sp.simplify(matdetfun[1][0])*((sp.simplify(matdetfun[0][1])*sp.simplify(matdetfun[2][2]))-(sp.simplify(matdetfun[2][1])*sp.simplify(matdetfun[0][2]))))+
            (sp.simplify(matdetfun[2][0])*((sp.simplify(matdetfun[0][1])*sp.simplify(matdetfun[1][2]))-(sp.simplify(matdetfun[1][1])*sp.simplify(matdetfun[0][2])))))

#-------------------------
#Determinant of NxN Matrix

def det_nxn(mat, total=0):
    index = list(range(len(mat)))
    if len(mat) == 2 and len(mat[0]) == 2:
        val = det2(mat)
        return val
    for i in index:
        mat1 = copymat(mat) 
        mat1 = mat1[1:]
        height = len(mat1)
        for j in range(height):
            mat1[j] = mat1[j][0:i] + mat1[j][i+1:]

        sign = (-1) ** (i % 2) 
        sub_det = det_nxn(mat1)
        total += sign * mat[0][i] * sub_det
    return total

#-----------------
#Def Main Function

print_mat(finalmat)
detval=sp.simplify(det_nxn(finalmat))
print()
print("Determinant of the Wronskian")
print()
print(detval)
print()
print("After Substituting Zero in the Det Equation",detval.subs(x, 0))

if(detval.subs(x, 0)!=0):
  print()
  print("The Given equations are Linearly Independent")
else:
  print()
  print("Linearly Dependent")