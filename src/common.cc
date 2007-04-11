#include <common.h>

void fatal_error(char *s) {
	fprintf(stderr, "error: %s\n", s);
	exit(1);
}

void *malloc_or_die(size_t size) {
	if(size <= 0) fatal_error("size <= 0 in malloc_or_die");
	void *p = malloc(size);
	if(!p) fatal_error("out of memory");
	return p;
}

void *realloc_or_die(void *p, size_t size) {
	if(size <= 0) fatal_error("size <= 0 in realloc_or_die");
	p = realloc(p, size);
	if(!p) fatal_error("out of memory");
	return p;
}
