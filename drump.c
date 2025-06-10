#include <time.h>
#include <pthread.h>
#include <sys/types.h>
#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"
#define BRIDGE 310
#define DOWN 0
#define UP 250
#define KEYS 17

/* here are our X variables */
Display *dis;

int screen;
int deadlock[20];
unsigned char *x_buffer;
unsigned char *image3;
Window win;
GC gc;
XImage *x_image;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();


struct node { int x; int y; double z; double zv; int nc; int nb[4];};
struct drum { int x; int y; int rate; int radius; int open; double *image; int tid; double *wav; int note[10000];};


void disp (double *,int,int);

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}

int qlott(double *image, int xp, int yp, int r, int g, int b, int xs, int ys)
{
	int x,y;
	if (xs<1){ xs=1;}
	if (ys<1){ ys=1;}
	for (x=-xs/2;x<=xs/2;x++)
	{
		for (y=-ys/2;y<=ys/2;y++)
		{
			plotf(image,xp+x,yp+y,r,g,b);
		}
	}
}	

void init_n(struct node **n ,struct node **m,int len)
{
	int i;
	for (i=0;i<len;i++)
	{
        	*(n+i)=(struct node *)malloc(sizeof (struct node)); //nodes 
        	*(m+i)=(struct node *)malloc(sizeof (struct node)); //nodes 
		n[i]->nc=0;
		m[i]->nc=0;
		n[i]->zv=0;
		m[i]->zv=0;
	}

}

void strike_n(struct node **n, int len , double radius)
{
	int i;
	double amp;
	amp=800+rand()%400;
	for (i=0;i<len;i++)
	{
		n[i]->z=amp-(amp*sqrt((n[i]->y)*(n[i]->y))+((n[i]->x)*(n[i]->x))/radius);
	}

}

void do_nm(struct node **n ,struct node **m, int tid, int len, double mass, double damping, int open)
{
	int i,j;
	double xa,d;

	d=0.9998+damping ;

	for (i=0;i<len;i++)
	{
			double force;
		if (n[i]->nc==4)
		{
			force=-(4*n[i]->z)+(n[n[i]->nb[0]]->z+n[n[i]->nb[1]]->z+n[n[i]->nb[2]]->z+n[n[i]->nb[3]]->z);
			m[i]->zv=(n[i]->zv*d)+(force/mass);
			m[i]->z=n[i]->z+(m[i]->zv);
		}else if (n[i]->nc>0 ){
		//}else if (n[i]->nc>0 ){
			if (open)
			{
			force=-((double)n[i]->nc*n[i]->z);
			for (j=0;j<n[i]->nc;j++)
			{
				force+=n[n[i]->nb[j]]->z;
			}
			m[i]->zv=(n[i]->zv*d)+(force/mass);
			m[i]->z=n[i]->z+(m[i]->zv);
			}else{
			m[i]->z=0;}
		} 

	}

}

void do_n(struct node **n ,struct node **m, int tid,int len,double mass, double damping, int open)
{
	do_nm(n ,m, tid,len,mass,damping,open);
	do_nm(m ,n, tid,len,mass,damping,open);
}

void draw_n (double *image, struct node **n,double xxp,double yyp,int len, int open)
{
	int i;
	int xp,yp;
	int sz;
	sz=4;
	for (i=0;i<len;i++)
	{
		xp=(n[i]->x*sz)+xxp; yp=(n[i]->y*sz)+yyp;
		if (n[i]->nc<4)
		{
			//if (bongo){qlott(image,xp,yp,0,255,255,14,14);}
			if (open){qlott(image,xp,yp,0,255,255,sz,sz);}
			else{qlott(image,xp,yp,255,0,255,sz,sz);}
		}else{
			qlott(image,xp,yp,127+(127*sin(n[i]->z/520)),
					127+(127*sin(n[i]->z/500)),
					127+(127*sin(n[i]->z/480)),sz,sz);
		}
	}
	//xp=(n[strike]->x*10)+(X_SIZE/2); yp=(n[strike]->y*10)+(Y_SIZE/2);
	//qlott(image,xp,yp,0,0,255,20,20);
	//xp=(n[left]->x*3)+(X_SIZE/2); yp=(n[left]->y*3)+(Y_SIZE/2);
	//qlott(image,xp,yp,255,0,255,20,20);
	//xp=(n[right]->x*3)+(X_SIZE/2); yp=(n[right]->y*3)+(Y_SIZE/2);
	//qlott(image,xp,yp,0,255,255,20,20);

}



void* skin(void* arg) {
	struct drum *dr;
	struct node **n,**m;
	int i,j,radius,r2,len,x,y,left,right;
	int rate;
	double damping;

	// Get our pointer
	dr=arg;
    	printf("Created a new thread radius %d\n",dr->radius);
	//struct drum { int x; int y; int rate; int radius; int open; double *image; };

	radius=dr->radius;
	r2=radius*radius;
	len=10*r2;
	rate=48000;

	n = (struct node**)malloc(len * sizeof(struct node*)); 
	m = (struct node**)malloc(len * sizeof(struct node*)); 
	init_n(n,m,len);

	// Build the shape
        len=0;
        for (x=-radius;x<radius;x++)
        {
                for (y=-radius;y<radius;y++)
                {
                        //if ((y*y)+((double)x*(double)x)<r2  && rand()%10000>40) // split
                        if ((y*y)+((double)x*(double)x)<r2 ) // split
                        {
                        n[len]->x=x;
                        n[len]->y=y;
                        n[len]->zv=0;
                        m[len]->zv=0;
                        m[len]->z=0;
                        len++;
                        }
                }
        }

        for (i=0;i<len;i++)
        {
                int ux,uy,dx,dy,lx,ly,rx,ry;
                ux=n[i]->x; uy=n[i]->y+1;
                dx=n[i]->x; dy=n[i]->y-1;
                lx=n[i]->x-1; ly=n[i]->y;
                rx=n[i]->x+1; ry=n[i]->y;
                for (j=0;j<len;j++)
                {
                        if ((n[j]->x==ux && n[j]->y==uy) || (n[j]->x==dx && n[j]->y==dy) ||
                            (n[j]->x==lx && n[j]->y==ly) || (n[j]->x==rx && n[j]->y==ry))
                        {
                                n[i]->nb[n[i]->nc]=j;
                                m[i]->nb[n[i]->nc]=j;
                                n[i]->nc++;
                                m[i]->nc++;
                        }

                }

        }
	left=(len/2)+(radius/2);
	right=(len/2)-(radius/2);

	double mass;
	mass=70+((dr->tid/2)*(dr->tid/2)*10);
	damping=0;
	int nc;
	nc=0;
	int framec,counter;
	counter=0;

	while (1)
	{
		for ( framec=0;framec<dr->rate;framec++)
		{
			if (counter%(rate/8)==0)
			{
				if (dr->note[nc]!=-1){
					damping=(double)dr->note[nc]*0.00001;
					strike_n(n,len ,radius);
					printf ("STRIKE %d %d\n",nc, dr->note[nc]);
				}
				nc++;
			}
			n[0]->z=0;
			do_n(n,m,dr->tid,len,mass,damping,dr->open);
			dr->wav[counter*2]=n[left]->z;
			dr->wav[(counter*2)+1]=n[right]->z;
			counter++;
		}
//void draw_n (double *image, struct node **n,double xxp,double yyp,int len)
		draw_n (dr->image,n,dr->x,dr->y,len,dr->open);
		deadlock[dr->tid]=0;
		while(deadlock[dr->tid]==0){ usleep (1000);}
	}
    	return NULL;
 }




int main(int argc,char *argv[])
{
	double *image2;
	double *wavs[20];
	
	int a,b,i,j,loop,frame,x,y;
	FILE *list;
	struct drum *dr[20];
	short *wav;
	long addup;
	int rate,dur,total,frames,drums;

	pthread_t thread1[20];

	rate=48000; dur=600; frames=30; drums=10;

        image2=(double *)malloc(sizeof (double)*X_SIZE*Y_SIZE*3); // disp buffer
        image3=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        wav=(short *)malloc(sizeof (short)*2*rate*dur); // wav buffer

	init_x();

    	// Creating a new thread.
	for (i=0;i<drums;i++)
	{
        	wavs[i]=(double *)malloc(sizeof (double)*2*rate*dur); // disp buffer
		dr[i] = (struct drum*)malloc(sizeof(struct drum)); 
		dr[i]->image=image2;
		dr[i]->x=(380*(i/2))+100;
		dr[i]->y=(Y_SIZE/4)+(i%2*Y_SIZE/2);
		dr[i]->radius=25+(4*(i/2));
		dr[i]->tid=i;
		dr[i]->open=i%2;
		dr[i]->wav=wavs[i];
		dr[i]->rate=rate/frames;
		deadlock[i]=1;
		for (j=0;j<10000;j++){ dr[i]->note[j]=-1;}
	}



  	if ((list= fopen("drump.lst", "r")) == NULL) {
    		fprintf(stderr, "can't open %s\n", "tune.lst");
    	return 0;}



	int pc;
	char line[300];
	pc=0;
    	while (fgets(line,300,list) != NULL) {
	    	int k;
	    	for (k=0;k<drums;k++)
	    	{
		    if (line[k]>47 && line[k]<58){ dr[k]->note[pc]=(int)line[k]-48;}
		    line[k]=0;
	    	}
		printf ("%s",line);
	    	printf ("skin 0 %d\n",dr[0]->note[pc]);
	    	printf ("skin 1 %d\n",dr[1]->note[pc]);
		pc++;
   	}

	for (i=0;i<drums;i++)
	{
		pthread_create(&thread1[i], NULL, skin, dr[i]);
	}

        printf("I found %d notes\n", pc);

	int go;
	frame=0; go=1;
	clearf(image2,0);

	while (go)
	{
		int tot;
		tot=0;
		for (i=0;i<drums;i++){ tot+=deadlock[i];}
		if (tot==0)
		{
			disp(image2,frame,1);
			frame++;
			for (i=0;i<drums;i++){ deadlock[i]=1;}
			total+=rate/frames;
			printf ("Doing Frame %d\n",frame);
			if (frame>=pc*frames/8){go=0;}
		}else 
		{
		       	usleep (1000);
		}
	}

	printf ("here \n");


	for (i=0;i<pc*rate/4;i++)
	{
		int p;
		p=0;
		for (j=0;j<drums;j++)
		{ 
			p+=dr[j]->wav[i]*10;
		}
		if (p>32767){p=32767;}
		if (p<-32767){p=-32767;}
		wav[i]=p;
	}

	save_wav(wav,"twang.wav",rate, 2, pc*rate/4 );

	char junk[30];
	scanf("%c",junk);
	close_x();
	exit(0);
}	

void disp (double *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*X_SIZE*3;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=255-image2[xpoint];
			x_buffer[X_POINT+1]=255-image2[xpoint+1];
			x_buffer[X_POINT]=255-image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/pl%05d.jpg",fram);
	if (ab){
		for (x=0;x<3*X_SIZE*Y_SIZE;x++){ image3[x]=255-image2[x];}
		jayit(image3,X_SIZE, Y_SIZE, input);}
	free (input);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

