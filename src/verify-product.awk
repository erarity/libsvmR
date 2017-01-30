#!/usr/bin/gawk -f
BEGIN{
    FS=" *|:"
}

/Skip list/{
    getline;
    for(i=1; i<=NF; i++){
	skips[$i] = $i;
    }
}
/Dataset/{
    getline;
    row = 0
    while(NF > 0){
	for(i=1; i<=NF; i+=2){
	    x[row][$i] = $(i+1)
	}
	row++;
	getline;
    }

    product = 0
    for(i=0; i<row; i++){
	for(j=0; j<row; j++){
	    for(k in x[i]){
		if(skips[k]==""){
		    product += x[i][k] * x[j][k]
		}
	    }
	}
    }
    print "Verified product " product;
}

/Sum/{
    print;
}

