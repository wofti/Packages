/* WT: if we set 
#define EXTERN extern
we can compile without -fcommon, while
#define EXTERN 
should result in the old behavior */
#define EXTERN extern
