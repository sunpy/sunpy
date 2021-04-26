Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/focal64"

  # Setup machine
  config.vm.provision "shell", privileged: true, inline: <<-SHELL
    # Set debconf interface
    export DEBIAN_FRONTEND=noninteractive
    # Update package lists
    apt-get -qq update > /dev/null
    # Install build tools (gcc, python3, pip3, etc)
    apt-get -qq install gcc python3-dev python3-pip python3-venv python3-wheel > /dev/null
  SHELL

  # Setup sunpy dev environment
  config.vm.provision "shell", privileged: false, inline: <<-SHELL
    cd /vagrant
    # Create virtual environment inside ~/.venv folder
    python3 -m venv ~/.venv
    # Install dependencies inside the virtual environment
    source ~/.venv/bin/activate
    pip3 install -q --upgrade pip
    python3 --version && pip3 list
    pip3 install -e ".[dev]"
  SHELL

  # Setup login shell
  config.vm.provision "shell", privileged: false, inline: <<-SHELL
    # Disable Ubuntu welcome message
    touch ~/.hushlogin
    # Autologin inside the virtual environment
    echo "source ~/.venv/bin/activate" >> /home/vagrant/.bashrc
    # Auto cd into shared folder
    echo "cd /vagrant" >> /home/vagrant/.bashrc
    # Print useful versions on login
    echo "python3 --version" >> /home/vagrant/.bashrc
  SHELL
end
