Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/bionic64"

  # Setup machine
  config.vm.provision "shell", privileged: true, inline: <<-SHELL
    apt install build-essential -y
  SHELL

  # Install miniconda
  config.vm.provision "shell", privileged: true, inline: <<-SHELL
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p /home/vagrant/miniconda
    rm miniconda.sh
    chown -R vagrant:vagrant /home/vagrant/miniconda
  SHELL
  config.vm.provision "shell", privileged: false, inline: <<-SHELL
    bash -i -c " \
      /home/vagrant/miniconda/bin/conda init bash && \
      source ~/.bashrc && \
      conda update -n base -c defaults conda \
    "
  SHELL

  # Setup sunpy dev environment
  config.vm.provision "shell", privileged: false, inline: <<-SHELL
    cd /vagrant
    bash -i -c " \
      conda create -y -c=conda-forge -n=sunpy-dev pip && \
      conda activate sunpy-dev && \
      pip install -e .[dev] \
    "
  SHELL

  # CD into shared folder and activate conda env on login
  config.vm.provision "shell", privileged: false, inline: <<-SHELL
    echo "cd /vagrant" >> /home/vagrant/.bashrc
    echo "conda activate sunpy-dev" >> /home/vagrant/.bashrc
    echo "conda --version && python --version" >> /home/vagrant/.bashrc
  SHELL
end
