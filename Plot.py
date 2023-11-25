import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle

# Import instance
with open('__Pics/Sol_T1_5_0.txt') as file:
    # First line contains the truck info 
    L,W,n = map(int, file.readline().strip().split(" "))
    
    # One line per item
    columns = []
    for i in range(n):
        column = [int(number) for number in file.readline().strip().split(" ")]
        print(column)
        columns.append(column)

# Get names
names = [str(columns[i][4]) for i in range(n)]
for i in range(n):
    for j in range(5,len(columns[i])):
        names[i] += "-" + str(columns[i][j])

# Draw instance
fig, ax = plt.subplots()
ax.plot([0,W,W,0,0],[0,0,L,L,0],color="black")
print(columns)

for i in range(n):
    ax.add_patch(Rectangle((columns[i][1], columns[i][0]), columns[i][3], columns[i][2], edgecolor='k', linewidth=1, facecolor='y'))
    ax.text(columns[i][1] + 0.5*columns[i][3], columns[i][0] + 0.5*columns[i][2], names[i],
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=5, color='black')
plt.show()
