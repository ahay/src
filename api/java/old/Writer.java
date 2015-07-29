package rsf;

import java.io.*;
import java.nio.ByteOrder;
import edu.mines.jtk.io.ArrayOutputStream;

public class Writer{

    public static void writeHeader(Header header,String fullPath){
            try{
                    PrintWriter pw = null;
                    
                    pw = new PrintWriter(new FileOutputStream(header.name));
                    pw.println("Warning:  Java programs were used!  This file may not be reproducible.");
                    pw.println();
                    pw.println("in=\""+fullPath+"\"");
                    pw.println("data_format="+"\""+header.format+"\"");
                    pw.println("esize="+header.esize);
                    for(int i = 0, j=1; i < Header.MAX_DIMS; ++i, ++j){
                        if (header.n[i] > 1){
                                pw.println("n"+j+"="+header.n[i]);
                                pw.println("d"+j+"="+header.deltas[i]);
                                pw.println("o"+j+"="+header.origins[i]);
                                if (header.labels[i] != null){
                                        pw.println("label"+j+"=\""+header.labels[i]+"\"");
                                } 
                                if (header.units[i] != null){
                                        pw.println("unit"+j+"=\""+header.units[i]+"\"");
                                }
                        }
                    }
                    pw.flush();
                    //pw.close();
            } catch (Exception e){
                    System.err.println(e);
            }             
    }

    public static void writeRSF(Header header, float[] data, String filename){
            header.name = filename;
            String fullPath = System.getenv("DATAPATH")+ header.name+"@";
            writeHeader(header,fullPath);
            //System.err.printf("\n\n%c%c%c",'\014','\014','\004');
            try{
                    ArrayOutputStream aos = null;
              //      if (fullPath.equals("stdin")){
                //        aos = new ArrayOutputStream(ps,ByteOrder.nativeOrder());
                  //  } else {
                        aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());
                    //}
                    aos.writeFloats(data);
                    aos.close();
            } catch (Exception e){
                    System.err.println(e);
                    }
    }

    public static void writeRSF(Header header, float[][] data, String filename){
            header.name = filename;
            String fullPath = System.getenv("DATAPATH")+ header.name+"@";
            writeHeader(header,fullPath);
           // System.err.printf("\n\n%c%c%c",'\014','\014','\004');
            try{
                    ArrayOutputStream aos = null;
                   
                    aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());
   
                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.err.println(e);
                    }
    }
    
    public static void writeRSF(Header header, float[][][] data, String filename){
            header.name = filename;
            String fullPath = System.getenv("DATAPATH")+ header.name+"@";
            writeHeader(header,fullPath);
            //System.err.printf("\n\n%c%c%c",'\014','\014','\004');
            try{
                    ArrayOutputStream aos = null;
                    aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());
                    aos.writeFloats(data);
                    aos.close();
            } catch (Exception e){
                    System.err.println(e);
                    }
    }
}

