## Test EggLib call from Terminal

# All markers

/Users/vitorpavinato/anaconda/bin/python summstats.py input-file=input_egglib_1.txt output-file=output_egglib_1.txt LSS=He,Dj,WCst WSS=He,Dj,WCst,S,thetaW,D,Da,ZZ GSS=He,Dj,WCst,S,thetaW,D,Da,ZZ wspan=1 select=all

# Random markers

/Users/vitorpavinato/anaconda/bin/python summstats.py input-file=input_egglib_1.txt output-file=output_egglib_2.txt LSS=He,Dj,WCst WSS=He,Dj,WCst,S,thetaW,D,Da,ZZ GSS=He,Dj,WCst,S,thetaW,D,Da,ZZ wspan=1 select=rand select-num=8

# Frequency of markers

/Users/vitorpavinato/anaconda/bin/python summstats.py input-file=input_egglib_1.txt output-file=output_egglib_3.txt LSS=He,Dj,WCst WSS=He,Dj,WCst,S,thetaW,D,Da,ZZ GSS=He,Dj,WCst,S,thetaW,D,Da,ZZ wspan=1 select=freq select-freq=2



# Only select=Y markers - one marker each chromosome

/Users/vitorpavinato/anaconda/bin/python summstats.py input-file=input_egglib_2.txt output-file=output_egglib_4.txt LSS=He,Dj,WCst WSS=He,Dj,WCst,S,thetaW,D,Da,ZZ GSS=He,Dj,WCst,S,thetaW,D,Da,ZZ wspan=1 select=list