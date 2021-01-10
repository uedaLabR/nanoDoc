path_r="C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\indexM6Amod.txt"
path_w="C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\indexM6Amod_2.txt"

with open(path_w, mode='w') as fw:
    f = open(path_r)
    line = True # 1行を文字列として読み込む(改行文字も含まれる)
    cnt = 0
    index = 0
    chrb4 = 0
    while line:
        line = f.readline()
        if(len(line)<10):
            break
        line = line.rstrip('\n')
        chr = line.split(",")[1]
        if (chr != chrb4) or (cnt%4000==0):
            index  = index + 1
        fw.write(line +"," +str(index)+"\n")
        cnt = cnt+1
        chrb4 = chr

    f.close

