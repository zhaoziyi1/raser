import os
import sys

def main():
    args = sys.argv[1:]
    n_number = int(args[0])
    t_number = int(args[1])
    n=int(t_number/n_number)
    for i in range(n):
        intance_n  = "instance"+str(i)
        os.system("singularity instance start raser.simg "+intance_n)

    for i in range(n):
        print(i)
        e_number=n_number*(i+1)
        os.system("nohup singularity exec instance://instance%s ./run 0.2.2 %s %s %s &" %(i,e_number,t_number,n_number))

if __name__ == '__main__':
    main()