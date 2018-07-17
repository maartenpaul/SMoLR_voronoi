def get_area(x):

   y=list(map(get_polygon_area,x))
   
   return y
  
#from https://www.geeksforgeeks.orgarea-of-a-polygon-with-given-n-ordered-vertices
def get_polygon_area(x):
    X = x[:,0]
    Y = x[:,1]
    # Initialze area
    area = 0.0
    n = len(X)
    # Calculate value of shoelace formula
    j = n - 1
    for i in range(0,n):
        area += (X[j] + X[i]) * (Y[j] - Y[i])
        j = i   # j is previous vertex to i
     
 
    # Return absolute value
    return abs(area / 2.0)
