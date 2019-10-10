""" 
Función para derivar en la dirección y tomando en cuenta la esfericidad de la tierra

INPUTS
array_2D: numpy array de dos dimensiones
dy: paso espacial en dirección y

OUTPUTS
dUy: array de dos dimensiones con las derivadas

""" 

#%%

# Funcion derivada segun y tomando en cuenta la esfericidad de la tierra
def derivy(array_2D,dy):
    # array_2D es el la matrix, dy el paso vertical
    dUy=array_2D*0
    a=len(array_2D[:,0])
    b=len(array_2D[0,:])
    
    i=0
    while i<a:
        j=0
        while j<b:
            if i==0:
                dUy[i,j]=(array_2D[i+1,j]-array_2D[i,j])/dy
            elif i==a-1:
                dUy[i,j]=(array_2D[i,j]-array_2D[i-1,j])/dy
            else:
                dUy[i,j]=(array_2D[i+1,j]-array_2D[i-1,j])/dy/2
            j=j+1
        i=i+1
    return dUy
