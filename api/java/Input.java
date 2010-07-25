package rsf;


public class Input extends RSFFile {

    public Input(String filename){
        super(filename,true);
    }

    public void read(float[] array){
        m8r.sf_floatread(array,(long)array.length,this.file);
    }
    
    public void read(float[][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][][][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][][][][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][][][][][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
    public void read(float[][][][][][][][][] array){
        for(int i = 0; i < array.length; ++i){
            read(array[i]);
        }
    } 
    
}
