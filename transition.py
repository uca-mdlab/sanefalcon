ls -1 | while read line; do NAME=$(echo $line | cut -d '.' -f1-3,7) ; mv $line $NAME ;done

