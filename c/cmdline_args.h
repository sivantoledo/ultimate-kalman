#ifndef CMDLINE_ARGS_H
#define CMDLINE_ARGS_H

void parse_args       (int argc, char *argv[]);
int  get_int_param    (const char* key, int*   value, int         default_value);
int  get_string_param (const char* key, char** value, const char* default_value);
int  get_boolean_param(const char* key, int*   value, int         default_value);
void check_unused_args(void);

#endif /* CMDLINE_ARGS_H */
