# fedora-live-m8r.ks
#
# Fedora Live Spin with m8r and the a minimal version of the XFCE Desktop
#
# Usage:
# 0. Have a machine with Fedora 15 x86_64
# 1. Install dependencies. As root:
#    yum -y install livecd-tools system-config-kickstart spin-kickstarts
# 2. Set SELinux to permissive mode. As root:
#    setenforce 0
# 3. Make the LiveCD:
#    livecd-creator -f m8rFedoraLive -c fedora-live-m8r.ks

%include /usr/share/spin-kickstarts/fedora-live-base.ks
%include /usr/share/spin-kickstarts/fedora-live-minimization.ks

repo --name=Madagascar --includepkgs=madagascar --baseurl=http://ahay.org/yum/repo/fedora/15/x86_64

%packages

# Graphics
evince
-evince-dvi
-evince-djvu
-gimp

# development
-geany
vim-enhanced
emacs

# Internet
firefox
-remmina

# More Desktop stuff
NetworkManager-vpnc
NetworkManager-openvpn
NetworkManager-gnome
NetworkManager-pptp
desktop-backgrounds-compat
-gparted
-setroubleshoot
xdg-user-dirs-gtk

# Command line
ntfs-3g
-powertop
wget

# xfce packages
@xfce-desktop
Terminal
gtk-xfce-engine
-ristretto
hal-storage-addon
-thunar-volman
-thunar-media-tags-plugin
gigolo
xfce4-battery-plugin
xfce4-cellmodem-plugin
xfce4-clipman-plugin
xfce4-cpugraph-plugin
xfce4-datetime-plugin
xfce4-diskperf-plugin
xfce4-fsguard-plugin
xfce4-genmon-plugin
xfce4-mount-plugin
xfce4-netload-plugin
xfce4-places-plugin
xfce4-power-manager
xfce4-quicklauncher-plugin
xfce4-screenshooter-plugin
xfce4-sensors-plugin
xfce4-smartbookmark-plugin
xfce4-systemload-plugin
xfce4-taskmanager
xfce4-time-out-plugin
xfce4-timer-plugin
xfce4-verve-plugin
xfce4-xfswitch-plugin
xfce4-xkb-plugin
-xfwm4-themes

-aspell-*
-man-pages-*
-gimp-help
-desktop-backgrounds-basic
-autofs
-nss_db
-acpid

# drop some system-config things
-system-config-boot
-system-config-language
-system-config-lvm
-system-config-network
-system-config-rootpassword
-system-config-services
-policycoreutils-gui

-iok
-sendmail
-jwhois
-finger
-smolt-firstboot
-abrt
-abrt-*
-orage
-rsync
-cdparanoia-libs
-cracklib*
-cvs
-alsa-*
-flac
-cyrus-*
-gnome-bluetooth-libs
-m17n-*
-kacst-*
-anthy
-kasumi

# Replace official Fedora logos with generics because this is an unofficial
# remix:
-fedora-logos
-fedora-bookmarks
-fedora-release* 
generic-logos
generic-release

# m8r dependencies
binutils
gcc
glibc-headers
scons
gcc-c++
gcc-gfortran
numpy
python
swig
libgomp
blas
blas-devel
atlas
atlas-devel
scipy
libXaw-devel

madagascar

%end

%post
# xfce configuration

# create /etc/sysconfig/desktop (needed for installation)

cat > /etc/sysconfig/desktop <<EOF
PREFERRED=/usr/bin/startxfce4
EOF

mkdir -p /root/.config/xfce4

cat > /root/.config/xfce4/helpers.rc <<EOF
MailReader=sylpheed-claws
EOF

cat >> /etc/rc.d/init.d/livesys << EOF

mkdir -p /home/liveuser/.config/xfce4

cat > /home/liveuser/.config/xfce4/helpers.rc << FOE
MailReader=sylpheed-claws
FOE

# disable screensaver locking
cat >> /home/liveuser/.xscreensaver << FOE
lock:	False
FOE
# set up timed auto-login for after 60 seconds
cat >> /etc/gdm/custom.conf << FOE
[daemon]
TimedLoginEnable=true
TimedLogin=liveuser
TimedLoginDelay=60
FOE

# Show harddisk install on the desktop
sed -i -e 's/NoDisplay=true/NoDisplay=false/' /usr/share/applications/liveinst.desktop
mkdir /home/liveuser/Desktop
cp /usr/share/applications/liveinst.desktop /home/liveuser/Desktop

# Madagascar configuration
echo 'export DATAPATH=/var/tmp/' >> /home/liveuser/.bashrc

# this goes at the end after all other changes. 
chown -R liveuser:liveuser /home/liveuser
restorecon -R /home/liveuser

EOF

%end
