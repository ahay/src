% replaces certain strings for the Madagascar manual



% \\texttt{
% x and y
% x or y


replace('xx and yy' ,'\\texttt{xx} and \\texttt{yy}','*.tex')
replace('x and y'   ,'\\texttt{x} and \\texttt{y}'  ,'*.tex')
replace('xx or yy'  ,'\\texttt{xx} or \\texttt{yy}' ,'*.tex')
replace('x or y'    ,'\\texttt{x} or \\texttt{y}'   ,'*.tex')
replace('size of xx','size of \\texttt{xx}'         ,'*.tex')
replace('size of x' ,'size of \\texttt{x}'          ,'*.tex')
replace('size of yy','size of \\texttt{yy}'         ,'*.tex')
replace('size of y' ,'size of \\texttt{y}'          ,'*.tex')

replace('whether adj is','whether \\texttt{adj} is','*.tex')
replace('whether add is','whether \\texttt{add} is','*.tex')


replace('_init}'  ,'\\_init}'  ,'*.tex')
replace('_define}','\\_define}','*.tex')
replace('_close}' ,'\\_close}' ,'*.tex')
replace('_solve}' ,'\\_solve}' ,'*.tex')
replace('\\\\_}'  ,'\\_}'      ,'*.tex')

%% Types
replace('type const float*'  ,'type \\texttt{const float*}'    ,'*.tex')
replace('type const float'   ,'type \\texttt{const float}'     ,'*.tex')
replace('type const int*'    ,'type \\texttt{const int*}'      ,'*.tex')
replace('type const int'     ,'type \\texttt{const int}'       ,'*.tex')
replace('type const char*'   ,'type \\texttt{const char*}'     ,'*.tex')
replace('type const char'    ,'type \\texttt{const char}'      ,'*.tex')


replace('type vc3d','type \\texttt{vc3d}','*.tex')
replace('type pt3d','type \\texttt{pt3d}','*.tex')


replace('type sf_operator2'   ,'type \\texttt{sf\\_operator2}'   ,'*.tex')
replace('type sf_operator'    ,'type \\texttt{sf\\_operator}'    ,'*.tex')
replace('type sf_stack'       ,'type \\texttt{sf\\_stack}'       ,'*.tex')
replace('type sf_complex***'  ,'type \\texttt{sf\\_complex***}'  ,'*.tex')
replace('type sf_complex**'   ,'type \\texttt{sf\\_complex**}'   ,'*.tex')
replace('type sf_complex*'    ,'type \\texttt{sf\\_complex*}'    ,'*.tex')
replace('type sf_complex'     ,'type \\texttt{sf\\_complex}'     ,'*.tex')
replace('type sf_file'        ,'type \\texttt{sf\\_file}'        ,'*.tex')
replace('type sf_axis'        ,'type \\texttt{sf\\_axis}'        ,'*.tex')
replace('type sf_solverstep'  ,'type \\texttt{sf\\_solverstep}'  ,'*.tex')
replace('type sf_csolverstep' ,'type \\texttt{sf\\_csolverstep}' ,'*.tex')
replace('type sf_coperator'   ,'type \\texttt{sf\\_coperator}'   ,'*.tex')
replace('type sf_interpolator','type \\texttt{sf\\_interpolator}','*.tex')
replace('type sf_eno'         ,'type \\texttt{sf\\_eno}'         ,'*.tex')
replace('type sf_eno2'        ,'type \\texttt{sf\\_eno2}'        ,'*.tex')
replace('type sf_filter'      ,'type \\texttt{sf\\_filter}'      ,'*.tex')
replace('type sf_list'        ,'type \\texttt{sf\\_list}'        ,'*.tex')

replace('type const sf_axis'    ,'type \\texttt{const sf\\_axis}'    ,'*.tex')
replace('type const sf_bands'   ,'type \\texttt{const sf\\_bands}'   ,'*.tex')
replace('type const sf_complex*','type \\texttt{const sf\\_complex*}','*.tex')


replace('type float***','type \\texttt{float***}','*.tex')
replace('type float**' ,'type \\texttt{float**}' ,'*.tex')
replace('type float*'  ,'type \\texttt{float*}'  ,'*.tex')
replace('type float'   ,'type \\texttt{float}'   ,'*.tex')

replace('type double***','type \\texttt{double***}','*.tex')
replace('type double**' ,'type \\texttt{double**}' ,'*.tex')
replace('type double*'  ,'type \\texttt{double*}'  ,'*.tex')
replace('type double'   ,'type \\texttt{double}'   ,'*.tex')

replace('type void*','type \\texttt{void*}','*.tex')
replace('type void' ,'type \\texttt{void}' ,'*.tex')


replace('type int***','type \\texttt{int***}','*.tex')
replace('type int**' ,'type \\texttt{int**}' ,'*.tex')
replace('type int*'  ,'type \\texttt{int*}'  ,'*.tex')
replace('type int'   ,'type \\texttt{int}'   ,'*.tex')


replace('type off_t*','type \\texttt{off\\_t*}','*.tex')
replace('type off_t' ,'type \\texttt{off\\_t}' ,'*.tex')
replace('type const off_t*','type \\texttt{const off\\_t*}','*.tex')
replace('type const off_t', 'type \\texttt{const off\\_t}' ,'*.tex')


replace('type size_t', 'type \\texttt{size\\_t}' ,'*.tex')

replace('type uint','type \\texttt{uint}','*.tex')
replace('type bool','type \\texttt{bool}','*.tex')
replace('type kiss_fft_cpx','type \\texttt{kiss\\_fft\\_cpx}','*.tex')



% replace('type sf_','type \\texttt{sf\\_}','*.tex')


%% 
% replace('\\item[\\texttt{','\\item[','*.tex')
replace('}]',']','*.tex')


% Mad_lib To do
% ===========
% '. Must be of (the) XXX type .' -->' (XXX).'
% '. Must be of (the) type XXX.' -->' (XXX).'
% datatypes
% 
% 
% . It is of type XXX. --> (XXX).
% Check Call and Definitions
% is computed.
replace('    ...', '   ...', '*.tex');

typs    = {'bool', 'char', 'double', 'FILE', 'float', 'int', ...
           'off\\_t', 'pt3d', 'pt2d', 'pt3d', ...
           'size\\_t', 'uchar', 'vc1d', 'vc2d', 'vc3d', 'void'};
tic
for tp = length(typs):length(typs)
    ast = ['****'];
    for i = 0:4
        dtyp = [typs{tp}, ast];
        fprintf('%s\n',dtyp);
        
        orig = ['. Must be of type \\texttt{',dtyp,'}.'];
%         orig = ['. Must be of the type \\texttt{',dtyp,'}.'];
%         orig = ['. Must be of the \\texttt{',dtyp,'} type.'];
%         orig = ['. Must be of \\texttt{',dtyp,'} type.'];
        new  = [' (\\texttt{', dtyp, '}).'];
%         fprintf('%s --> %s\n',orig,new);
        replace( orig, new, '*.tex');

        orig = ['. Must be of type \\texttt{const ',dtyp,'}.'];
%         orig = ['. Must be of the type \\texttt{const ',dtyp,'}.'];
%         orig = ['. Must be of the \\texttt{const ',dtyp,'} type.'];
%         orig = ['. Must be of \\texttt{const ',dtyp,'} type.'];
        new  = [' (\\texttt{const ', dtyp, '}).'];
%         fprintf('%s, %s\n',orig,new);
        replace( orig, new, '*.tex');
        
        ast = ast(1:end-1);
    end
end

sf_typs = {'axis', 'bands', 'complex', 'coperator', 'csolverstep', ...
           'list', 'operator2','operator','stack', 'solverstep',...
           'datatype', 'eno', 'eno2', 'file', 'filter', 'interpolator'};
        
for tp = 1:length(sf_typs)
    ast = ['****'];
    for i = 0:4
        dtyp = [sf_typs{tp}, ast];
        fprintf('%s\n',dtyp);

%         orig = ['. Must be of type \\texttt{sf\\_',dtyp,'}.'];
%         orig = ['. Must be of the type \\texttt{sf\\_',dtyp,'}.'];
%         orig = ['. Must be of the \\texttt{sf\\_',dtyp,'} type.'];
        orig = ['. Must be of \\texttt{sf\\_',dtyp,'} type.'];
        new  = [' (\\texttt{sf\\_', dtyp, '}).'];
%         fprintf('%s, %s\n',orig,new);
        replace( orig, new, '*.tex');

%         orig = ['. Must be of type \\texttt{const sf\\_',dtyp,'}.'];
%         orig = ['. Must be of the type \\texttt{const sf\\_',dtyp,'}.'];
%         orig = ['. Must be of the \\texttt{const sf\\_',dtyp,'} type.'];
        orig = ['. Must be of \\texttt{const sf\\_',dtyp,'} type.'];
        new  = [' (\\texttt{const sf\\_', dtyp, '}).'];
        replace( ['const ',orig], new, '*.tex');
        
        ast = ast(1:end-1);
    end
end
toc

%%

replace('true or flase','true or false','*.tex')


replace('\\subsection{','\\subsection{\\texttt{','*.tex')
replace('\\subsection{\\texttt{','\\subsection{{','*.tex')
replace('bool value','boolean value','*.tex')


replace('\\begin{list}','\\begin{itemize}','*.tex')
replace('\\end{list}','\\end{itemize}','*.tex')

       
   
replace('$ is computed.','$.','*.tex')
replace('{\\tt }{\\quad}[','{\\tt }{\\quad}[\\tt ','*.tex')

replace('   \\item[','   \\item[','*.tex')
