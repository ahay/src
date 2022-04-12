安装curvelops全过程
处理流程在/geo/yangliu/apefsnsep/curvemod

1）安装python3
https://baijiahao.baidu.com/s?id=1720756753507481024&wfr=spider&for=pc
第一步：安装相关依赖包和编译环境
$yum -y install libffi-devel
$yum -y install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel tk-devel gdbm-devel db4-devel libpcap-devel xz-devel gcc

（注意：这一步很重要，如果不安装相关依赖包，在使用pip安装python包时会报找不到SSL错误！）

第二步：下载python3.8.0安装包

$wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tar.xz

第三步：解压安装包并创建安装目录

$xz -d Python-3.8.0.tar.xz

$tar -xvf Python-3.8.0.tar

$mkdir /usr/local/python3.8.0

第四步：编译安装

$cd Python-3.8.0

$./configure --with-ssl --prefix=/usr/local/python3.8.0

(注意：prefix后面配置第三步中创建的路径，且等号两边不能有空格，不然会报错；如果用这个路径则需要root用户，如果用/home/yangliu/geotools/Python-3.8.0则不用root用户)

$make && make install

第五步：创建python3.8软链接

$ln -s /usr/local/python3.8.0/bin/python3.8 /usr/bin/python3
注意修改路径
$ln -s /usr/local/python3.8.0/bin/pip3.8 /usr/bin/pip3
注意修改路径
$pip3 install --upgrade pip（升级pip3）

第六步：修改python2.7.5软链接（这一步可有可无）

$mv /usr/bin/python /usr/bin/python2

第七步：验证，使用python3进入python3.8.0命令行

$python3

2）安装fftw-2.1.5
$ wget http://www.fftw.org/fftw-2.1.5.tar.gz
$ tar -xvf fftw-2.1.5.tar.gz
$ cd fftw-2.1.5/
$ make clean
$  ./configure --with-pic --prefix=/home/yangliu/geotools/fftw-2.1.5 --with-gcc=/usr/bin/gcc
$ make
$ make install
$ cd ..

3）安装CurveLab
下载安装包从http://www.curvelet.org/（需要注册）或者直接使用下载文件
$ tar -xvf CurveLab-2.1.3.tar.gz
$ cd CurveLab-2.1.3
$ make clean
$ emacs makefile.opt
修改其中FFTW路径为FFTW_DIR = 	/home/yangliu/geotools/fftw-2.1.5
保存退出makefile.opt文件
$ make lib
$ make test （出现时间等计算信息）
$ cd ..

4）修改用户shell文件
$ emacs ~/.bashrc
export FFTW=/home/yangliu/geotools/fftw-2.1.5
export FDCT=/home/yangliu/geotools/CurveLab-2.1.3
保存退出

5）安装各种python库函数
$ pip3 install pylops
$ pip3 install pybind11
$ pip3 install PyWavelets

6）安装curvelops
（1）一种方法是按照官网直接下载安装，但是可能遇到github无法访问的情况
python3 -m pip install git+https://github.com/cako/curvelops@main
（2）另一种方法是找到安装包curvelops-main，再安装
$ cd curvelops-main/
如果python3是用root用户安装则需要以下三步
$ python3 setup.py build
$ su
$ python3 setup.py install
如果python3是普通用户安装则需要以下一步
$ python3 setup.py install
测试是否安装完成
$ python3
>>>from curvelops import FDCT2D
>>>
则表示安装正常

7）测试处理流程
$ cd ~/geo/yangliu/apefsnsep/curvemod
$ scons
此时会报错，是因为部分程序用到了curvelet_denoise2d.py的ct_denoise.rsf结果，不用管
$ python3 curvelet_denoise2d.py
此时生成ct_denoise.rsf
$ scons lock

