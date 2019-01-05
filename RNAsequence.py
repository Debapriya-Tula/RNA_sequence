# Mappings
main_dict = {"A":('U',),"G":('C','U'),
             "U":('A','G'),"C":('G',)}


#ACGU - ACGU format of storing in the 2d list
AU = [[0.0, 0.0, 0.0, -0.9],
      [0.0, 0.0, -2.2, 0.0],
      [0.0, -2.1, 0.0, -0.6],
      [-1.1, 0.0, -1.4, 0.0],
     ]

CG = [[0.0, 0.0, 0.0, -2.1],
      [0.0, 0.0, -3.3, 0.0],
      [0.0, -2.4, 0.0, -1.4],
      [-2.1, 0.0, -2.1, 0.0],
     ]

GC = [[0.0, 0.0, 0.0, -0.9],
      [0.0, 0.0, -2.2, 0.0],
      [0.0, -2.1, 0.0, -0.6],
      [-1.1, 0.0, -1.4, 0.0],
     ]

GU = [[0.0, 0.0, 0.0, -1.3],
      [0.0, 0.0, -2.5, 0.0],
      [0.0, -2.1, 0.0, -0.5],
      [-1.4, 0.0, 1.3, 0.0],
     ]

UA = [[0.0, 0.0, 0.0, -1.3],
      [0.0, 0.0, -2.4, 0.0],
      [0.0, -2.1, 0.0, -1.0],
      [-0.9, 0.0, -1.3, 0.0],
     ]

UG = [[0.0, 0.0, 0.0, -1.0],
      [0.0, 0.0, -1.5, 0.0],
      [0.0, -1.4, 0.0, 0.3],
      [-0.6, 0.0, -0.5, 0.0],
     ]

#*******************************************************************************



#Printing the matrix
def print_matrix(arr):
    l1 = len(arr)
    l2 = len(arr[0])

    for i in range(l1):
        for j in range(l2):
            print(arr[i][j],end = "     ")
        print()

'''
# For bifurcation
def find_max(c,i,j):
    maxm = 0
    t = i+1
    for k in range(i+1,j):
        if (c[i][k] + c[k+1][j]) > maxm:
            maxm = c[i][k] + c[k+1][j]
            t = k
    return maxm,t
'''



#The algorithm that gives the matrix to finding the optimal structure
def find_stem_loop(arr,l1):

    #Initialisation of the matrices with the required values of gaps
    c = [[0 for j in range(l1)] for i in range(l1)]

    m = 0; n = 1
    while n < l1:
        i = m; j = n
        while j < l1:
            # Choosing the max among the three values of pairings
            x = c[i+1][j-1] + (1 if arr[j] in main_dict[arr[i]] else 0)
            y = c[i+1][j]
            z = c[i][j-1]

            #w,_ = find_max(c,i,j) # For bifurcation
            c[i][j] = max(x,y,z)
            i += 1
            j += 1
        n += 1

    return c




#take input from file
f1 = open("data.txt","r")
r = f1.readlines()
f1.close()

r = r[1]
l1 = len(r)-1
arr = str(r[0:l1])
print("The sequence taken from file:",r)

# The matrix for traceback
result_matrix = find_stem_loop(arr,l1)

print("The resultant matrix:")
print_matrix(result_matrix)


#*******************************************************************************



'''
class Bifur:
    def __init__(self,i,j,k):
        self.data1 = (i,k)
        self.data2 = (k+1,j)
        self.b_left = None
        self.b_right = None
'''

# For storing the paths in a ternary-tree
class Node:
    # constructor to create tree node
    def __init__(self,i,j):
        self.data = (i,j)
        self.left = None
        self.mid = None
        self.right = None


    # finding all paths that could be possible optimal structures
    def reconstruct(self,arr,i,j,c):
        if c[i][j] == 0:
            return
        else:
            # temp variable to store the value derived from the diagonally opposite value
            temp = c[i+1][j-1] + (1 if arr[j] in main_dict[arr[i]] else 0)

            #bifur,k = find_max(c,i,j)
            t = max(temp,c[i+1][j],c[i][j-1])

            if t == c[i][j-1]:
                self.left = Node(i,j-1)
                self.left.reconstruct(arr,i,j-1,c)

            if t == temp:
                self.mid = Node(i+1,j-1)
                self.mid.reconstruct(arr,i+1,j-1,c)

            if t == c[i+1][j]:
                self.right = Node(i+1,j)
                self.right.reconstruct(arr,i+1,j,c)

        return self



# Function for printing the paths
def printPaths(root):
    # list to store path
    path = []
    printPathsRec(root, path, 0)



# Helper function to print path from root to leaf in binary tree
def printPathsRec(root, path, pathLen):
    if root is None:
        return

    if(len(path) > pathLen):
        path[pathLen] = root.data
    else:
        path.append(root.data)

    pathLen = pathLen + 1

    if root.left is None and root.mid is None and root.right is None:
        # leaf node then print the list
        printArray(path, pathLen)

    else:
        # try for left,mid and right subtree
        printPathsRec(root.left, path, pathLen)
        printPathsRec(root.mid, path, pathLen)
        printPathsRec(root.right, path, pathLen)



#Storing all possible paths in an list named 'store'
store = []
def printArray(ints, len):
    store.append(ints[0 : len])


# Calling the reconstruct function from the top-right corner element in result_matrix as root
root = Node(0,l1-1)
root = root.reconstruct(arr,0,l1-1,result_matrix)

printPaths(root)
print("\nThe paths:")
print_matrix(store)


#*******************************************************************************



# the find_loc function assigns the appropriate energy value acc. to the tables listed above
def find_loc(en,a,b):
    dic = {"A":0,"C":1,"G":2,"U":3}
    return en[dic[a]][dic[b]]

# access the right table for energy considerations
def cal_energy(arr,i,j):
    t1 = arr[i+1]
    t2 = arr[j-1]

    energy = 0
    if arr[i]=="A" and arr[j]=="U":
        energy = find_loc(AU,t1,t2)

    elif arr[i]=="C" and arr[j]=="G":
        energy = find_loc(CG,t1,t2)

    elif arr[i]=="G":
        if arr[j]=="C":
            energy = find_loc(GC,t1,t2)
        elif arr[j]=="U":
            energy = find_loc(GU,t1,t2)

    elif arr[i]=="U":
        if arr[j]=="A":
            energy = find_loc(UA,t1,t2)
        elif arr[j]=="U":
            energy = find_loc(UG,t1,t2)

    return energy



#Print  the paths in ".()" format
l2 = len(store)
final = [["." for j in range(l1)] for i in range(l2)]

cal = [0 for i in range(l2)]
for i in range(l2):
    for j in range(len(store[i])-1):
        t = store[i][j]
        t1 = store[i][j+1]

                                                # This is the condition for a stack(t1[0]...t[1]-1)
        if arr[t[0]] in main_dict[arr[t[1]]] and (t1[0]==t[0]+1 and t1[1]==t[1]-1) and abs(t[0]-t[1])>2:
            final[i][t[1]] = ")"
            final[i][t[0]] = "("
            # call the cal_energy fn. to assign a value only if there is a stack formed
            cal[i] += cal_energy(arr,t[0],t[1])

    # Evaluate for the end-point in the result_matrix for a possible structure
    j += 1
    t = store[i][j]
    if arr[t[0]] in main_dict[arr[t[1]]] and abs(t[0]-t[1])>2:
        final[i][t[1]] = ")"
        final[i][t[0]] = "("



# For sorting the list1 on the basis of list2 in reverse order
def sort_list(list1, list2):
    zipped_pairs = zip(list2, list1)
    z = [x for _, x in sorted(zipped_pairs)]
    return z



# Write output to an output file
f2 = open("write.txt","w")
final = sort_list(final,cal)
final = ["".join(i) for i in final]
f2.writelines("Info about the output\n")

for i in range(l2):
    f2.writelines(final[i]+"\n")


#*******************************************************************************
