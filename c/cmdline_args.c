/*
 * cmdline_args.c
 *
 * Simple but effective parsing of key=value command-line arguments.
 *
 * Code generated using Grok 3.
 *
 * (C) Sivan Toledo, 2025
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "cmdline_args.h"

#define MAX_ARGS      16
#define MAX_KEY_LEN   32
#define MAX_VALUE_LEN 32

typedef struct {
  char key[MAX_KEY_LEN];
  char value[MAX_VALUE_LEN];
  int used;
} keyvalue_t;

static keyvalue_t args[MAX_ARGS];
static int arg_count = 0;

void parse_args(int argc, char *argv[]) {
  arg_count = 0;
  for (int i = 1; i < argc && arg_count < MAX_ARGS; i++) {
    char *equals = strchr(argv[i], '=');
    if (equals && equals > argv[i] && equals < argv[i] + strlen(argv[i]) - 1) {
      size_t key_len = equals - argv[i];
      if (key_len >= MAX_KEY_LEN) {
        fprintf(stderr, "Key too long: %.*s\n", (int) key_len, argv[i]);
        exit(1);
      }
      strncpy(args[arg_count].key, argv[i], key_len);
      args[arg_count].key[key_len] = '\0';
      const char *value = equals + 1;
      if (strlen(value) >= MAX_VALUE_LEN) {
        fprintf(stderr, "Value too long for key: %s\n", args[arg_count].key);
        exit(1);
      }
      strcpy(args[arg_count].value, value);
      args[arg_count].used = 0;
      arg_count++;
    } else {
      fprintf(stderr, "Invalid argument format: %s (expected key=value)\n", argv[i]);
      exit(1);
    }
  }
}

int get_int_param(const char *key, int *value, int default_value) {
  for (int i = 0; i < arg_count; i++) {
    if (strcmp(args[i].key, key) == 0) {
      char *endptr;
      long val = strtol(args[i].value, &endptr, 10);
      if (*endptr == '\0') {
        *value = (int) val;
        args[i].used = 1;
        return 1;
      } else {
        fprintf(stderr, "Parameter %s must be an integer, got: %s\n", key, args[i].value);
        exit(1);
      }
    }
  }
  *value = default_value;
  return 0;
}

int get_string_param(const char *key, char **value, const char *default_value) {
  for (int i = 0; i < arg_count; i++) {
    if (strcmp(args[i].key, key) == 0) {
      *value = args[i].value;
      args[i].used = 1;
      return 1;
    }
  }
  *value = (char*) default_value;
  return 0;
}

int get_boolean_param(const char *key, int *value, int default_value) {
  for (int i = 0; i < arg_count; i++) {
    if (strcmp(args[i].key, key) == 0) {
      // Convert value to lowercase for case-insensitive comparison
      char lower_value[MAX_VALUE_LEN];
      strncpy(lower_value, args[i].value, MAX_VALUE_LEN - 1);
      lower_value[MAX_VALUE_LEN - 1] = '\0';
      for (char *p = lower_value; *p; p++) {
        *p = tolower((unsigned char) *p);
      }

      if (strcmp(lower_value, "true") == 0 || strcmp(lower_value, "1") == 0) {
        *value = 1;
        args[i].used = 1;
        return 1;
      } else if (strcmp(lower_value, "false") == 0 || strcmp(lower_value, "0") == 0) {
        *value = 0;
        args[i].used = 1;
        return 1;
      } else {
        fprintf(stderr, "Parameter %s must be a boolean (true/false/1/0), got: %s\n", key, args[i].value);
        exit(1);
      }
    }
  }
  *value = default_value;
  return 0;
}

void check_unused_args(void) {
  for (int i = 0; i < arg_count; i++) {
    if (!args[i].used) {
      fprintf(stderr, "Unexpected or unrecognized argument: %s=%s\n", args[i].key, args[i].value);
      exit(1);
    }
  }
}
