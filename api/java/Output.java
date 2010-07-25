package rsf;

public class Output extends RSFFile{

    public Output(String filename){
        super(filename,false);
    }
    
    public void write(float[] data){
            try {
                m8r.sf_floatwrite(data,(long)data.length,this.file);
            } catch (Exception e){
                System.err.println("FILE IO-------- ERROR ------\n");
                System.err.println(e);
            }            
    }

    public void write(float[][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public void write(float[][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public void write(float[][][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public void write(float[][][][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public void write( float[][][][][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public void write( float[][][][][][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public  void write( float[][][][][][][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
    
    public void write(float[][][][][][][][][] data){
            try{
                    for(int i = 0; i < data.length; ++i){
                        write(data[i]);
                    }
            } catch (Exception e){
                    System.err.println(e);
            }
    }
}

