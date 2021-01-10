import matplotlib.pyplot as plt

p = "C:\\Users\\hirok\\Desktop\\covid19\\index.txt"
f = open(p)
lines = f.readlines()
x = []
y = []
for line in lines:
    data = line.split(",")
    start = int(data[2])
    end = int(data[3])
    x.append(start)
    y.append(end)
f.close()

plt.scatter(x,y,s=1)
plt.show()