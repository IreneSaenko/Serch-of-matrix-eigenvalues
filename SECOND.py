import random
import copy

def normir(vect1):#нормирует вектор
    norma=0
    vect=vect1
    for i in range(N):
         norma+=vect[i]*vect[i]
    for i in range(N):
        vect[i]=vect[i]/pow(norma,0.5);
    return vect

def print_matrix(lst):#печатает матрицу
    print('matrix:\n')
    for i in range(N):
        print(lst[i],'\n')

def a_lambdae(lst,sigma):#вычитает sigma из главной диагонали
    for i in range(N):
        lst[i][i]=lst[i][i]-sigma
    return(lst)

def gauss(a,b):#решает систему методом Гаусса
 
 for k in range(N-1):#номера строк, вычетаемых из всех остальных
   
    for i in range(k+1,N):
         if abs(a[k][k])<1e-5:#если всё плохо, меняем ведущий элемент
             mx=k
             for z in range(k,N):
                if abs(a[z][k])>abs(a[mx][k]):
                    mx=z
             for n in range(N):
                    t=a[k][n]
                    a[k][n]=a[mx][n]
                    a[mx][n]=t
             t=b[k]
             b[k]=b[mx]
             b[mx]=t
    
         t = a[i][k]/a[k][k]
         for j in range(N):
              a[i][j] -= a[k][j]*t
         b[i] -= b[k]*t

 y=['*' for i in range(N)]
 i=N-1
 while i>=0:
        k=b[i]
        for j in range(i+1,N):
            k-=a[i][j]*y[j];
        y[i]=k/a[i][i]
        i-=1;
 return(y)


def multimatrix(vect,matrix):#умножает матрицу на вектор
            res=[]
            for i in range(N):
                t=0
                for j in range(N):
                    t+=vect[j]*matrix[i][j]
                res.append(t)
            return res

def norma(vect):#находит норму вектора
    norma=0
    for i in range(N):
         norma+=vect[i]*vect[i]
    return pow(norma,0.5);

def max(vect):#находит максимальную компоненту вектора(против переполнения)
    mx=0
    for i in range(N):
            if vect[i]> vect[mx]:
                mx=i
    return mx

        
#-------------------Начало алгоритма, поиск приближений степенным методом------
 
eps=0.3

f=open('matrix1.txt','r')
N=int(f.read(1))#считываем размерность

l=[]#список, хранящий матрицу
mainmatrix=[]
for line in f:
        mainmatrix.append([float(x) for x in line.split()])
        l.append([float(x) for x in line.split()])
l.remove([])
mainmatrix.remove([])

print_matrix(l)#выводим матрицу




x=[1,0,0]#ищем начальные приближения при помощи степенного метода

z1=0
alf1=norma(x)+3*eps
alf2=norma(x)
while (abs(alf2-alf1))>eps:
    alf1=alf2
    z1=x
    y=multimatrix(x,l)
    alf2=norma(y)
    x=normir(y)

lam1=alf2#начальное приближение старшего по модулю собственного числа
z2=multimatrix(x,l)
j=max(x)
lam2=(z2[j]-lam1*x[j])/(x[j]-z1[j])

print('Начальные приближения, найденные степенным методом:\n')
print('lambda1=',lam1,'\n')
print('lambda2=',lam2,'\n')



#--------------------------Начало метода обратных итераций------------------


eps=1e-7

print('Точность:',eps,'\n')

x=[]#задание случайного вектора для начального приближения 
for i in range(N):
    x.append(random.randint(1,10))
x=normir(x)


l1=a_lambdae(l,lam1)#нашли матрицу A-lambda*E для старшего lambda


z1=0
z2=3*eps
while abs(z2-z1)>eps:
    s=copy.deepcopy(l1)
    t=x
    y=gauss(s,x)
    x=normir(y)
    j=max(t)
    z1=z2
    z2=x[j]/t[j]



print('Собственный вектор, соответствующий lambda1:\n')
print(x)


x=[1,1,1]
x=normir(x)

l2=a_lambdae(mainmatrix,lam2)#нашли матрицу A-lambda*E для старшего lambda


z1=0
z2=3*eps
while abs(z2-z1)>eps:
    s=copy.deepcopy(l2)
    t=x
    y=gauss(s,x)
    x=normir(y)
    j=max(t)
    z1=z2
    z2=x[j]/t[j]

print('Собственный вектор, соответствующий lambda2:\n')
print(x)



