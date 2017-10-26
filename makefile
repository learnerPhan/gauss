test: readlinebyline.c
	gcc -g -std=c99 -o test readlinebyline.c

clean:
	rm ./test
