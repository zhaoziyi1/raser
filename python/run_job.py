import os
import sys

def main():
    args = sys.argv[1:]
    n_number = int(args[0])
    t_number = int(args[1])
    n=int(t_number/n_number)
    for i in range(n):
        e_number=n_number*(i+1)
        os.system("nohup ./run 0.2.2 %s %s %s &" %(e_number,t_number,n_number))

if __name__ == '__main__':
    main()