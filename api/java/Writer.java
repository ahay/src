package rsf;

import java.io.*;
import java.nio.ByteOrder;
import edu.mines.jtk.io.ArrayOutputStream;

public class Writer{

    public static PrintStream writeHeader(Header header,String fullPath){
            PrintStream ps = null;
            try{
                    PrintWriter pw = null;
                    
                    if (header.name.equals("out")){
                        ps = new PrintStream(System.out);
                        pw = new PrintWriter(ps);
                    } else {
                        pw = new PrintWriter(new FileOutputStream(header.name));
                    }
                    String path = System.getenv("DATAPATH");
                    pw.println("java programs warning, this file may not be reproducibile");
                    pw.println();
                    pw.println("in=\""+fullPath+"\"");
                    pw.println("data_format="+"\""+header.format+"\"");
                    pw.println("esize="+header.esize);
                    for(int i = 0, j=1; i < header.ndims; ++i, ++j){
                        pw.println("n"+j+"="+header.n[i]);
                        pw.println("d"+j+"="+header.deltas[i]);
                        pw.println("o"+j+"="+header.origins[i]);
                        pw.println("label"+j+"=\""+header.labels[i]+"\"");
                        pw.println("unit"+j+"=\""+header.units[i]+"\"");
                    }
                    pw.flush();
                    //pw.close();
            } catch (Exception e){

                    System.out.println(e);
            }
            return ps;
                 
    }

    public static void writeRSF(Header header, float[] data, String filename){
            header.name = filename;
            String fullPath = "";
            if (filename.equals("out")){
                fullPath = "stdin";
            } else {
                fullPath = System.getenv("DATAPATH")+ header.name+"@";
            }
            PrintStream ps = writeHeader(header,fullPath);
            System.out.printf("\n\n%c%c%c",'\014','\014','\004');
            try{
                    ArrayOutputStream aos = null;
                    if (fullPath.equals("stdin")){
                        aos = new ArrayOutputStream(ps,ByteOrder.nativeOrder());
                    } else {
                        aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());
                    }
                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.out.println(e);
                    }
            ps.close();
    }

    public static void writeRSF(Header header, float[][] data, String filename){
            header.name = filename;
            String fullPath = "";
            if (filename.equals("out")){
                fullPath = "stdin";
            } else {
                fullPath = System.getenv("DATAPATH")+ header.name+"@";
            }
            PrintStream ps = writeHeader(header,fullPath);
            System.out.printf("\n\n%c%c%c",'\014','\014','\004');
            try{
                    ArrayOutputStream aos = null;
                    if (fullPath.equals("stdin")){
                        aos = new ArrayOutputStream(ps,ByteOrder.nativeOrder());
                    } else {
                        aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());
                    }
                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.out.println(e);
                    }
            ps.close();
    }
    public static void writeRSF(Header header, float[][][] data, String filename){
            header.name = filename;
            String fullPath = "";
            if (filename.equals("out")){
                fullPath = "stdin";
            } else {
                fullPath = System.getenv("DATAPATH")+ header.name+"@";
            }
            PrintStream ps = writeHeader(header,fullPath);
            System.out.printf("\n\n%c%c%c",'\014','\014','\004');
            try{
                    ArrayOutputStream aos = null;
                    if (fullPath.equals("stdin")){
                        aos = new ArrayOutputStream(ps,ByteOrder.nativeOrder());
                    } else {
                        aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());
                    }
                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.out.println(e);
                    }
            ps.close();


    }


}

