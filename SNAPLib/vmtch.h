#pragma once
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <limits.h>
#include <inttypes.h>
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <fcntl.h>
#include <search.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/param.h>
#include <sys/user.h>
#include <sys/sysctl.h>

#ifdef __linux__
#include <sys/ioctl.h>
#include <sys/mount.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

int64_t total_pages = 0;
int64_t total_pages_in_core = 0;
void vmcheck(char *fpath);
void vmtouch(char *path, bool check);

#ifdef __cplusplus
}
#endif
