#!/bin/bash

#--------------------------------------------------------

usage()
{
    cat << EOF
This script uses vagrant to start a virtual machine and 
builds percolator with converters packages by executing 
a builder file.

usage: $0 
                  [[-h]] [[-a]]
                  [[-b branch]]|[[-s source_directory]]
                  [[-r release_directory]]
                  -p ubuntu|centos|fedora|nw32|nw64

If no branch and source_directory is provided, the source
code from which the sourcecode is checked out from will be used.
Make sure that Vagrant and VirtualBox are up to date.
  -h     prints this help page
  -a     keeps the vagrant box alive (i.e. do not call vagrant destroy)
EOF
}

#--------------------------------------------------------
#input options management
#--------------------------------------------------------
#Default values:

current_path=$PWD;
script_dir=$(dirname ${BASH_SOURCE});
cd ${script_dir};
builder_adr="../builders/"
# boxes, might be overridden later on
vagbox_name="box-cutter/fedora23"
vagbox_url=""

while getopts “hab:s:r:p:” OPTION; do
    case $OPTION in
        h)  usage;exit 1;;
        a)  alive="1";;
        b)  branch=$OPTARG;;
        s)  src=$OPTARG;;
        r)  release=$OPTARG;;
        p)  case $OPTARG in
               	ubuntu)
                    post="ubuntu64"
                    vagbox_name="ubuntu/trusty64"
                    vagbox_url=""
                    package_ext="deb"
                    ;;
                fedora)	
                    post="fedora64"
                    package_ext="rpm"
                    ;;
                centos) 
                    post="centos64"
                    vagbox_name="bento/centos-7.2"
                    #vagbox_name="bento/centos-6.7"
                    vagbox_url=""
                    package_ext="rpm"
                    ;;
                nw32) 
                    post="nativew32"
                    batfile=true
                    vagbox_name="win10vs15"
                    vagbox_url="~/VagrantWin7/win10vs15.box"
                    package_ext="exe"
                    ;;
                nw64) 
                    post="nativew64"
                    batfile=true
                    vagbox_name="win10vs15"
                    vagbox_url="~/VagrantWin7/win10vs15.box"
                    package_ext="exe"
                      ;;
                *)
                    if [[ $OPTARG == *,* ]]; then
                      arr=$(echo $OPTARG | tr "," "\n");
                      multi_platform="1";
                    else
                      echo "Platform $OPTARG is undefined."
                      exit 1
                    fi;;
            esac;;
        \?)  echo "Invalid option: -${OPTARG}" >&2;;
    esac
done
if [[ -z $batfile ]]; then
  builder="${post}_build.sh";
else
  builder="${post}_build.bat";
fi

######
if [[ ! -z $multi_platform ]]; then
  [[ ! -z $branch ]] && arg=$arg"-b $branch "; 
  [[ ! -z $src ]] && arg=$arg"-s $src "; 
  [[ ! -z $release ]] && arg=$arg"-r $release "; 
  for x in $arr; do
    call_arg=$arg"-p $x";
    echo "Recusive call: $0 $call_arg"
    bash $0 $call_arg
  done
  exit 0;
fi
if [[ -z $post ]]; then
    usage
    echo "Please select one or more platforms with -p option."
    exit 1
fi
######

echo "------------------------------------------------------------------------";
echo "About to start up a build note using the builder script ${builder}";
echo "------------------------------------------------------------------------";


tmp_dir="$(mktemp -d --tmpdir tmp_${post}_XXXX)"
mkdir -p ${tmp_dir}/src/maracluster
if [[ -z $src ]]; then
  if [[ -z $branch ]]; then
    echo "Copying source code from ${script_dir} to ${tmp_dir}/src/maracluster/"
    cp -R ${script_dir}/../../* ${tmp_dir}/src/maracluster/
  else
    echo "Cloning source code using the branch ${branch}"
    git clone --branch ${branch} https://github.com/statisticalbiotechnology/maracluster.git ${tmp_dir}/src/maracluster
  fi
else
  echo "Copying source code from user specified path ${src}"
  cp -R ${src}/* ${tmp_dir}/src/maracluster
fi
######
if [[ -z $release ]]; then
  mkdir ${tmp_dir}/${post}_release
  release="${tmp_dir}/${post}_release"
fi

echo "Executing build procedure using:"
echo "tmp_dir=${tmp_dir}" 
echo "post=${post}" 
echo "tmp_dir=${tmp_dir}" 


#--------------------------------------------------------
#########################################################
#--------------------------------------------------------
#copy builder:
cp ${builder_adr}${builder} ${tmp_dir};
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# making the Vagrantfile:
cd ${tmp_dir};
touch Vagrantfile;

if [[ -z $batfile ]]; then

vagbox_url_line=""
if [[ -n ${vagbox_url} ]]; then
  vagbox_url_line="config.vm.box_url = \"${vagbox_url}\""
fi

#-----------------Vagrantfile content---------------
cat <<EOF > Vagrantfile
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "${vagbox_name}"
  ${vagbox_url_line}
  config.ssh.insert_key = false
  config.vm.boot_timeout = 600
  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "4096", "--cpus", "4"]
    # vb.gui = true # turn on for trouble shooting, e.g. if boot times out repeatedly
  end
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    shell.args = "-s /vagrant/src -r /vagrant/"
  end
end
EOF
#-----------------end of Vagrantfile content--------
#  config.vm.provision :shell, :inline => "su vagrant -c 'bash /vagrant/${builder} /vagrant/src /vagrant/build_${post}'"
else
cat <<EOF > Vagrantfile
Vagrant.configure("2") do |config|

  # Configure base box parameters
  config.vm.box = "${vagbox_name}"
  config.vm.box_url = "${vagbox_url}"
  config.vm.guest = :windows
  config.windows.halt_timeout = 30
  config.vm.boot_timeout = 1200

  # Port forward WinRM and RDP
  config.vm.communicator = "winrm"
  
  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--memory", "8192", "--cpus", "4"]
    # vb.gui = true # turn on for trouble shooting, e.g. if boot times out repeatedly
  end
  
  config.vm.provision :shell do |shell|
    shell.path = "${tmp_dir}/${builder}"
    shell.args = '-s "C:\vagrant\src" -b "C:\vagrant\build" -r "C:\vagrant"'
  end
end
EOF
fi
#---------------------------------------------------------------------------------------
vagrant up

#---------------------------------------------------------------------------------------
# release:

echo "Copying ready made packages from ${tmp_dir} to ${release}" 

mkdir -p ${release};
cp -v ${tmp_dir}/mara*.${package_ext} ${release};

#---------------------------------------------------------------------------------------

if [[ $? -eq 0 ]] && [[ -z ${alive} ]]; then
  vagrant destroy -f
else
  echo "-a option set or encountered error: keeping the VM alive, remember to close and delete the VM manually."
fi

#---------------------------------------------------------------------------------------
