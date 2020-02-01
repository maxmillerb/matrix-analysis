#ALGORITIMO PARA ANÁLISE DE PÓRTICO PLANO (PROC DESLOCAMENTOS)

#Este programa utiliza o numpy, pacote de Python que suporta operações com vetores e matrizes
import numpy as np

#LEITURA DOS DADOS

#Número de nós
Nnos=(int)(input('Digite o número de nós: '))
print('Há',Nnos,'nós')
#Número de elementos
Nelem=(int)(input('Digite o número de elementos: '))
print('Há', Nelem,'elementos')

#Coordenadas dos nós
X=np.zeros((1,Nnos))
Y=np.zeros((1,Nnos))
for i in range(Nnos):
    print('Coordenadas do nó',i+1,':')
    X[0,i]=(float)(input('Digite a coordenada x do nó: '))
    Y[0,i]=(float)(input('Digite a coordenada y do nó: '))
for i in range(0,Nnos):
    print('Nó',i+1,': (',X[0,i],', ',Y[0,i],')')

#Condições de contorno
#deslocamentos prescritos = 1
#deslocamentos desconhecidos = 0
cod=np.zeros((Nnos*3))
for i in range(Nnos):
    print('Condições de contorno do nó',i+1,':')
    cod[3*i]=(int)(input('O nó é livre na horizontal? 0=Sim ou 1=Não? '))
    cod[3*i+1]=(int)(input('O nó é livre na vertical? 0=Sim ou 1=Não? '))
    cod[3*i+2]=(int)(input('O nó é livre na rotação? 0=Sim ou 1=Não? '))

#Forças nodais atuantes
Fn=np.zeros((Nnos*3))
for i in range(Nnos):
    print('Forças concentradas atuantes no nó',i+1,':')
    Fn[3*i]=(float)(input('Qual o valor da força horizontal aplicado no nó? '))
    Fn[3*i+1]=(float)(input('Qual o valor da força vertical aplicado no nó? '))
    Fn[3*i+2]=(float)(input('Qual o valor do momento aplicado no nó? '))

#Conectividade dos elementos
con=np.zeros((Nelem,2))
for i in range(0,Nelem):
    print('O elemento',i+1,'conecta quais nós?')
    con[i,0]=(int)(input('Nó inicial: '))
    con[i,1]=(int)(input('Nó final: '))
print(con)

#Propriedades dos elementos
EA=np.zeros((1,Nelem))
EI=np.zeros((1,Nelem))
for i in range(0,Nelem):
    print('Propriedades do elemento',i+1,': ')
    EA[0,i]=(float)(input('Digite o valor de E*A do elemento: '))
    EI[0,i]=(float)(input('Digite o valor de E*I do elemento: '))
print('Vetor com os valores de EA de cada elemento: ', EA)
print('Vetor com os valores de EI de cada elemento: ', EI)

#Carregamentos distribuídos ao longo dos elementos (sistema de coordenadas globais)
qx=np.zeros((Nelem))
qy=np.zeros((Nelem))
for i in range(0,Nelem):
    print('Componentes da carga distribuída no elemento',i+1,': ')
    qx[i]=(float)(input('Digite a componente horizontal da carga distribuída no elemento: '))
    qy[i]=(float)(input('Digite a componente vertical da carga distribuída no elemento: '))

#Comprimentos dos elementos
L=np.zeros((Nelem))
con=con.astype(int)
for i in range(0,Nelem):
    j=con[i,0]
    k=con[i,1]
    xj=X[0,j-1]
    yj=Y[0,j-1]
    xk=X[0,k-1]
    yk=Y[0,k-1]
    L[i]=((xk-xj)**2+(yk-yj)**2)**0.5
    print('Comprimento do elemento',i+1,': ',L[i])

#Montagem da matriz de rigidez dos elementos no sistema local
rei=np.zeros((Nelem,6,6))
for i in range(0,Nelem):
    rei[i,0,0]=EA[0,i]/L[i]
    rei[i,0,3]=EA[0,i]/L[i]
    rei[i,3,3]=EA[0,i]/L[i]
    rei[i,3,0]=EA[0,i]/L[i]
    rei[i,1,1]=12*EI[0,i]/(L[i])**3
    rei[i,1,4]=12*EI[0,i]/(L[i])**3
    rei[i,4,1]=12*EI[0,i]/(L[i])**3
    rei[i,4,4]=12*EI[0,i]/(L[i])**3
    rei[i,2,2]=4*EI[0,i]/L[i]
    rei[i,5,5]=4*EI[0,i]/L[i]
    rei[i,1,5]=6*EI[0,i]/(L[i])**2
    rei[i,5,1]=6*EI[0,i]/(L[i])**2
    rei[i,4,5]=6*EI[0,i]/(L[i])**2
    rei[i,5,4]=6*EI[0,i]/(L[i])**2
    rei[i,1,2]=-6*EI[0,i]/(L[i])**2
    rei[i,2,1]=-6*EI[0,i]/(L[i])**2
    rei[i,2,4]=-6*EI[0,i]/(L[i])**2
    rei[i,4,2]=-6*EI[0,i]/(L[i])**2
    rei[i,2,5]=-2*EI[0,i]/(L[i])
    rei[i,5,2]=-2*EI[0,i]/(L[i])
    print('Matriz de rigidez do elemento',i+1,' no sistema local: ')
    print(rei[i])

#Montagem da matriz de incidência cinemática do elemento
beta=np.zeros((Nelem,6,6))
cos=np.zeros((Nelem))
sen=np.zeros((Nelem))

for i in range(0,Nelem):
    j=con[i,0]
    k=con[i,1]
    xj=X[0,j-1]
    yj=Y[0,j-1]
    xk=X[0,k-1]
    yk=Y[0,k-1]
    cos[i]=(xk-xj)/L[i]
    sen[i]=(yk-yj)/L[i]
    beta[i,0,0]=-cos[i]
    beta[i,4,4]=-cos[i]
    beta[i,1,1]=cos[i]
    beta[i,3,3]=cos[i]
    beta[i,0,1]=-sen[i]
    beta[i,1,0]=-sen[i]
    beta[i,3,4]=sen[i]
    beta[i,4,3]=sen[i]
    beta[i,2,2]=-1
    beta[i,5,5]=1
    print('Cosseno do ângulo de inclinação do elemento',i+1,': ')
    print(cos[i])
    print('Seno do ângulo de inclinação do elemento',i+1,': ')
    print(sen[i])
    print('Matriz de incidência cinemática do elemento',i+1,': ')
    print(beta[i])

#Cálculo da matriz de rigidez dos elementos no sistema global
betaT=np.zeros((Nelem,6,6))
rgi=np.zeros((Nelem,6,6))
for i in range(0,Nelem):
    betaT[i]=np.transpose(beta[i])
for i in range(0,Nelem):
    rgi[i]=np.matmul(np.matmul(betaT[i],rei[i]),beta[i])
    print('Matriz de rigidez do elemento',i+1,' no sistema global: ')
    print(rgi[i])

#Montagem da matriz de rigidez global da estrutura
R=np.zeros((3*Nnos,3*Nnos))
for i in range(Nelem):
    for j in range(2):
        for k in range(3):
            glib1=3*con[i,j]-3+k
            for l in range(2):
                for m in range(3):
                    glib2=3*con[i,l]-3+m
                    R[glib1,glib2]=R[glib1,glib2]+rgi[i,3*j+k,3*l+m]
print('Matriz de rigidez da estrutura: ')
print(R)

#Carregamentos distribuídos dos elementos (sistema de coordenadas local)
qxe=np.zeros((Nelem))
qye=np.zeros((Nelem))
for i in range(Nelem):
    qxe[i]=qx[i]*cos[i]+qy[i]*sen[i]
    qye[i]=qy[i]*cos[i]-qx[i]*sen[i]
print('qxe: ')
print(qxe)
print('qye: ')
print(qye)

#Reações nodais a carregamentos distribuídos
Pne_ei=np.zeros((Nelem,6))
for i in range(Nelem):
    Pne_ei[i,0]=qxe[i]*L[i]/2
    Pne_ei[i,1]=qye[i]*L[i]/2
    Pne_ei[i,2]=-qye[i]*L[i]**2/12
    Pne_ei[i,3]=-qxe[i]*L[i]/2
    Pne_ei[i,4]=-qye[i]*L[i]/2
    Pne_ei[i,5]=-qye[i]*L[i]**2/12
print('Pne_ei: ')
print(Pne_ei)

#Ações equivalentes dos carregamentos distribuídos dos elementos (sistema de coordenadas global)
Pne_gi=np.zeros((Nelem,6))
for i in range(Nelem):
    Pne_gi[i]=np.matmul(betaT[i],Pne_ei[i])
print('Pne_gi: ')
print(Pne_gi)

#Forças equivalentes
Fne=np.zeros((Nnos*3))
for i in range(Nelem):
    glib=3*(con[i,0])-3
    Fne[glib]=Fne[glib]+Pne_gi[i,0]
    Fne[glib+1]=Fne[glib+1]+Pne_gi[i,1]
    Fne[glib+2]=Fne[glib+2]+Pne_gi[i,2]
    glib=3*(con[i,1])-3
    Fne[glib]=Fne[glib]+Pne_gi[i,3]
    Fne[glib+1]=Fne[glib+1]+Pne_gi[i,4]
    Fne[glib+2]=Fne[glib+2]+Pne_gi[i,5]
print('Fne: ')
print(Fne)

#União das forças nodais com as forças equivalentes
F=np.zeros((Nnos*3))
F=Fn+Fne
print('F:')
print(F)

#Método de penalty
Rp=np.zeros((3*Nnos,3*Nnos))
nmg=10**10
for i in range(3*Nnos):
    for j in range(3*Nnos):
        Rp[i,j]=R[i,j]
        if cod[i]==1:
            Rp[i,i]=nmg
print('Matriz de rigidez da estrutura após o método de penalty:')
print(Rp)

#Cálculo dos deslocamentos globais (resolução do sistema de equações)
invRp=np.zeros((3*Nnos,3*Nnos))
invRp=np.linalg.inv(Rp)
D=np.zeros((3*Nnos,3*Nnos))
D=np.matmul(invRp,F)
print('Deslocamentos globais:')
print(D)

#Extração dos deslocamentos dos elementos (sistema global de coordenadas)
Dgi=np.zeros((Nelem,6))
for i in range(Nelem):
    Dgi[i,0]=D[3*(con[i,0]-1)]
    Dgi[i,1]=D[3*(con[i,0]-1)+1]
    Dgi[i,2]=D[3*(con[i,0]-1)+2]
    Dgi[i,3]=D[3*(con[i,1]-1)]
    Dgi[i,4]=D[3*(con[i,1]-1)+1]
    Dgi[i,5]=D[3*(con[i,1]-1)+2]

#Deslocamentos dos elementos (sistema local de coordenadas)
Dei=np.zeros((Nelem,6))
for i in range(Nelem):
    Dei[i]=np.matmul(beta[i],Dgi[i])

#Cálculo das reações
Reac=np.zeros((3*Nnos))
Reac=np.matmul(R,D)-Fne
print('Reações')
print(Reac)
for i in range(Nnos):
    print('Reação horizontal (positivo para a direita) no nó',i+1,':')
    print(Reac[3*i])
    print('Reação vertical (positivo para cima) no nó',i+1,':')
    print(Reac[3*i+1])
    print('Reação ao giro (positivo no sentido anti-horário) no nó',i+1,':')
    print(Reac[3*i+2])

#Esforços solicitantes nodais no elemento
Pei=np.zeros((Nelem,6))
for i in range(Nelem):
    Pei[i]=np.matmul(rei[i],Dei[i])-Pne_ei[i]
print('Esforços solicitantes nodais:')
print(Pei)
for i in range(Nelem):
    print('Esforço normal no extremo esquerdo (positivo para a esquerda) da barra',i+1,':')
    print(Pei[i,0])
    print('Esforço cortante no extremo esquerdo (positivo para cima) da barra',i+1,':')
    print(Pei[i,1])
    print('Momento fletor no extremo esquerdo (positivo no sentido horário) da barra',i+1,':')
    print(Pei[i,2])
    print('Esforço normal no extremo direito (positivo para a direita) da barra',i+1,':')
    print(Pei[i,3])
    print('Esforço cortante no extremo direito (positivo para baixo) da barra',i+1,':')
    print(Pei[i,4])
    print('Momento fletor no extremo direito (positivo no sentido anti-horário) da barra',i+1,':')
    print(Pei[i,5])

print('Fim')
