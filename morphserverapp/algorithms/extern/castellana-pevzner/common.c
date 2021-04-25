#include "common.h"

FILE *logfile;
const char log_name[] = "lib/bio/extern/castellana-pevzner/k.log";


void init_log() {
    logfile = fopen(log_name,"w");
    assert(logfile != NULL); 
    fclose(logfile);
}

void write_to_log(char *message) {
    logfile = fopen(log_name, "a");
    assert(logfile != NULL);
    fprintf(logfile, message);
    fclose(logfile);
}

