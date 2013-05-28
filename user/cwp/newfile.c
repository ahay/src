#include <rsf.h>
#include <sys/stat.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
/*^*/

static const int tabsize=10;

struct sf_File {
    FILE *stream, *head;
    char *dataname, *buf, *headname;
    sf_simtab pars;
    XDR xdr;
    enum xdr_op op;
    sf_datatype type;
    sf_dataform form;
    bool pipe, rw, dryrun;
};

static sf_file *infiles = NULL;
static size_t nfile=0;

off_t sf_byte (sf_file file)
/*< Count the file data size (in bytes) >*/
{
    int st;
    off_t size;
    struct stat buf;

    if (0 == strcmp(file->dataname,"stdin")) return ((off_t) -1);

    if (NULL == file->dataname) {
        st = fstat(fileno(file->stream),&buf);
    } else {
        st = stat(file->dataname,&buf);
    }
    if (0 != st) sf_error ("%s: cannot find file size:",__FILE__);
    size = buf.st_size;

    return size;
}

/*------------------------------------------------------------*/
sf_file sf_tmpfile(char *format)
/*< Create an temporary (rw mode) file structure. Lives within the program >*/
{ 

    sf_file file;
    char *dataname=NULL;
    off_t len;

    file = (sf_file) sf_alloc(1,sizeof(*file));
    file->stream = sf_tempfile(&dataname,"w+b");
        len = strlen(dataname)+1;
        file->dataname = (char*) sf_alloc(len,sizeof(char));
        memcpy(file->dataname,dataname,len);

    if (NULL == file->stream)
        sf_error ("%s: Trouble reading data file %s:",__FILE__,dataname);

    file->buf = NULL;
    file->pars = sf_simtab_init (tabsize);
    file->head = NULL;

    sf_putstring(file,"in",file->dataname);

    if (NULL == infiles) {
        infiles = (sf_file *) sf_alloc(1,sizeof(sf_file));
        infiles[0] = NULL;
        nfile=1;
    }

    if (NULL != infiles[0] &&
        NULL != (format = sf_histstring(infiles[0],"data_format"))) {
        sf_setformat(file,format);
        free (format);
    } else {
        sf_setformat(file,"native_float");
    }

    file->headname = file->dataname;
    file->dataname = NULL;

    return file;
}

/*------------------------------------------------------------*/
void sf_filefresh(sf_file file)
/*< used for temporary file only to recover the dataname >*/
{

    extern int fseeko(FILE *stream, off_t offset, int whence);

    file->dataname = file->headname;

    if (0 > fseeko(file->stream,(off_t)0,SEEK_SET))
        sf_error ("%s: seek problem:",__FILE__);
}


void sf_filecopy(sf_file file, sf_file src, sf_datatype type)
/*< copy the content in src->stream to file->stream >*/
{
    off_t nleft, n;
    char buf[BUFSIZ];
    extern int fseeko(FILE *stream, off_t offset, int whence);

    n = sf_bytes(src);
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);

    for (nleft = BUFSIZ; n > 0; n -= nleft) {
        if (nleft > n) nleft=n;
        switch (type) {
            case SF_FLOAT:
                sf_floatread  ((float*) buf,
                               nleft/sizeof(float),src);
                break;
            case SF_COMPLEX:
                sf_complexread((sf_complex*) buf,
                               nleft/sizeof(sf_complex),src);
                break;
            default:
                sf_error("%s: unsupported type %d",__FILE__,type);
                break;
        }
        if (nleft != fwrite (buf,1,nleft,file->stream))
            sf_error("%s: writing error:",__FILE__);
    }
}

void sf_tmpfileclose (sf_file file)
/*< close a file and free allocated space >*/
{
    if (NULL == file) return;

    if (file->stream != stdin &&
        file->stream != stdout &&
        file->stream != NULL) {
        (void) fclose (file->stream);
        file->stream = NULL;
    }

    if (file->headname != NULL) {
        (void) unlink (file->headname);
        free(file->headname);
        file->headname = NULL;
    }

    if (NULL != file->pars) {
        sf_simtab_close (file->pars);
        file->pars = NULL;
    }

    if (NULL != file->buf) {
        free (file->buf);
        file->buf = NULL;
    }
}



