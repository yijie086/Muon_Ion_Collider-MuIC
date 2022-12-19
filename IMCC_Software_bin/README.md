# Tutorial for Install Muon Collider Software on macOS

[toc]

## 1. Framework on macOS

### Dependency: Docker Desktop

Docker Desktop is an easy-to-install application for your Mac, Linux, or Windows environment that enables you to build and share containerized applications and microservices.

It provides a simple interface that enables you to manage your containers, applications, and images directly from your machine without having to use the CLI to perform core actions.

#### Download

For Mac with intel chip, the download link is: https://desktop.docker.com/mac/main/amd64/Docker.dmg?utm_source=docker&utm_medium=webreferral&utm_campaign=docs-driven-download-mac-amd64

- macOS must be version 10.15 or newer. That is, Catalina, Big Sur, or Monterey. We recommend upgrading to the latest version of macOS.
- At least 4 GB of RAM.
- VirtualBox prior to version 4.3.30 must not be installed as it is not compatible with Docker Desktop.

#### Installation

1. Double-click `Docker.dmg` to open the installer, then drag the Docker icon to the Applications folder.
2. Double-click `Docker.app` in the **Applications** folder to start Docker.
3. The Docker menu (![whale menu](https://docs.docker.com/desktop/install/images/whale-x.svg)) displays the Docker Subscription Service Agreement window.

### Dependency: XQuartz - an X11 display server for Mac OS X

XQuartz **allows cross-platform applications using X11 for the GUI to run on macOS**, many of which are not specifically designed for macOS. This includes numerous scientific and academic software projects.

#### Download

The down link is: https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.2/XQuartz-2.8.2.dmg

- version: 2.8.2
- For macOS 10.9 or later

#### Installation

1. Double-click `.dmg` file to open the installer, then double-click `.pkg` file.

## 2. Docker Pull the Image

Open `terminal.app` and run:

```
docker pull infnpd/mucoll-ilc-framework:1.6-centos8
```

## 3. Set the X11 display server

Run the `XQuartz.app` , then in it's application menu you'll find a Preferences choice. In the Preferences-Security make this setting:

```
Authenticate connections --> OFF
Allow connections from network clients --> YES
```

Then, restart `Xquartz` (or your whole machine).

The first thing is determining your IP address, this magic shell command gets the IP address into an environment variable:

```
export IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
```

```
echo $IP
```

With `XQuartz` running ... execute these commands:

```
xhost +$IP
```

```
defaults write org.xquartz.X11 enable_iglx -bool true
```

And then, keep `XQuartz` running.

## 4. Run Docker Containers

In Mac `terminal.app`, run the following command:

```
docker run -it --rm -e "DISPLAY=$IP:0" -e "XAUTHORITY=/.Xauthority" -e "USER=$USER" -w "/home/$USER" -v "/tmp/.X11-unix:/tmp/.X11-unix" -v "$HOME/.Xauthority:/home/$USER/.Xauthority" --net=host infnpd/mucoll-ilc-framework:1.6-centos8
```

### Explanation:

- `-it`: run in interactive mode, enter the container to view the content
- `--rm`: delete after use, if you want to keep the container, don't use this
- `-e`: specify environment variables
  - `-e "DISPLAY=$IP:0" -e "XAUTHORITY=/.Xauthority"`: set up the X11 display server
- `-w`:  working directory inside the container
- `-v`: mount the specified path
  - `-v "/tmp/.X11-unix:/tmp/.X11-unix" -v "$HOME/.Xauthority:/home/$USER/.Xauthority"`: set up the X11 display server
- `--net=host`: tell Docker not to put the container network into an isolated namespace, i.e. not to containerize the network inside the container. At this point the container uses localhost's network, which has full access to the localhost interface

After Running the beyond command, we enter the docker containers:

```
[root@docker-desktop #]
```

The container has the Muon Collider software installed, and we can start testing directly after sourcing the environment setting:

```
source /opt/ilcsoft/init_ilcsoft.sh
```

### Test 1: the X11 display server

Run following command in containers:

```
glxinfo
```

If X11 is working well, it will return something like this:

```
144 GLX Visuals
    visual  x   bf lv rg d st  colorbuffer  sr ax dp st accumbuffer  ms  cav
  id dep cl sp  sz l  ci b ro  r  g  b  a F gb bf th cl  r  g  b  a ns b eat
----------------------------------------------------------------------------
0x022 24 tc  0  32  0 r  y .   8  8  8  8 .  .  0 16  8  0  0  0  0  0 0 None
0x141 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0  0  0  0  0  0  0  0 0 Slow
0x142 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0  0  0  0  0  0  0 16 1 Slow
0x143 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0 32  0  0  0  0  0  0 0 Slow
0x144 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0 32  0  0  0  0  0 16 1 Slow
0x145 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0  0  8  0  0  0  0  0 0 Slow
0x146 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0  0  8  0  0  0  0 16 1 Slow
0x147 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0 32  8  0  0  0  0  0 0 Slow
0x148 24 tc  0  32  0 r  . .   8  8  8  8 .  .  0 32  8  0  0  0  0 16 1 Slow
0x149 24 tc  0  32  0 r  y .   8  8  8  8 .  .  0  0  0  0  0  0  0  0 0 Slow

............... AND LOTS OF OTHER SCREEN OUTPUT ............... 

```

### Test 2: the Muon Collider Software

Move to the test code's directory:

```
cd /opt/MuonCutil/SoftCheck/
```

#### Simulation

```
ddsim --compactFile /opt/ilcsoft/muonc/detector-simulation/geometries/MuColl_v1/MuColl_v1.xml --steeringFile sim_steer.py > sim.out
```

#### Reconstruction

```
Marlin --InitDD4hep_mod4.DD4hepXMLFile=/opt/ilcsoft/muonc/detector-simulation/geometries/MuColl_v1/MuColl_v1.xml reco_steer.xml > reco.out
```

#### Event display

```
ced2go -d /opt/ilcsoft/muonc/detector-simulation/geometries/MuColl_v1/MuColl_v1.xml Output_REC.slcio
```

Wait for a while, you can see the output graph on your screen.

At first, it will show something like this:

<img src="images/image-20221103175707677.png" alt="image-20221103175707677" style="zoom:30%;" />

And wait for several minutes, it will show the hole figure:

<img src="images/image-20221103175840089.png" alt="image-20221103175840089" style="zoom:37%;" />

<img src="images/image-20221103175905542.png" alt="image-20221103175905542" style="zoom:40%;" />

In the containers, double click for picking. Press `<ENTER>`  for the next event.

Inactivity for a long time will result in interruption of graphical output, so restart it in the container.

## 5. Beam Induced Background for Muon Collider

This tutorial focus on running simulation locally, so will download the data base from `login.snowmass21.io`.

### Download data base from `login.snowmass21.io`

First, you have to require membership from `https://connect.snowmass21.io`, using official mail with `.edu`.

#### Sign up

You should input your profile like the following:

```
Unix Username --> <ANY NAME YOU LIKE>		#will become your username when you login via ssh
Name --> <YOUR NAME>
Phone --> <YOUR PHONE>
Institution --> <YOUR INSTITUTION>
Email --> <YOUR INSTITUTIONAL EMAIL>
```

#### SSH public key

1. In a terminal, type:

   ```
   ssh-keygen -t rsa
   ```

2. Hit `enter` for the default location, and optionally enter a password. This will generate two files: A private key file (typically `id_rsa`) and a key file (typically `id_rsa.pub`). The private key should **never** be shared, and Snowmass21 Connect will never ask you for your private key.

3. In order to see your SSH public key type:

   `cat ~/.ssh/id_rsa.pub`

4. Use your mouse to select everything that is printed to the screen, the format should look like:

   `ssh-rsa AAAAB3N....M7Q== yourusername@yourmachine`

5. Copy the selection to the clipboard.

6. Paste the contents of the clipboard in the corresponding box on your CI Connect Profile.

When your account is active, you'll see this:

<img src="images/image-20221104111433864.png" alt="image-20221104111433864" style="zoom:30%;" />

#### Login via ssh

In `terminal.app`:

```
ssh -XY <Your Unix Username>@login.snowmass21.io
```

If there is no error information when you log in, that means X11 is working well, and you can using graphics UI if you want to run the code remotely on `login.snowmass21.io`. But, I don't recommend if you want to modify the simulation file.

Now `exit`:

```
exit
```

#### Download database

```
scp -r <Your Unix Username>@login.snowmass21.io:/collab/project/snowmass21/data/muonc <A local file path in your computer to store the data base>
```

### A look at the BIB

In Mac `terminal.app`, run the following command to create a container with the data base:

```
docker run -it --rm -e "DISPLAY=$IP:0" -e "XAUTHORITY=/.Xauthority" -e "USER=$USER" -w "/home/$USER" -v "/tmp/.X11-unix:/tmp/.X11-unix" -v "$HOME/.Xauthority:/home/$USER/.Xauthority" -v "<the local file path in your computer to store the data base>:/data" --net=host infnpd/mucoll-ilc-framework:1.6-centos8
```

Check the X11:

```
glxinfo
```

Sourcing the environment setting:

```
source /opt/ilcsoft/init_ilcsoft.sh
```

Check out the configuration files and the ROOT macros:

```
git clone https://github.com/MuonColliderSoft/MuC-Tutorial.git
```

```
cd  MuC-Tutorial/tutorial/2-BIB && git checkout v2.0
```

#### Beam-induced background (BIB) characterization

##### Fill the histograms

```
root -l make_mcPartPlots_BIB.C+
```

##### BIB particle entry point into the detector

```
h01_entryPoint->Draw("COLZ")
```

<img src="images/image-20221104113700021.png" alt="image-20221104113700021" style="zoom:40%;" />

#####  Z coordinates of entry point for different particles species

```
h02_entryZ_ph->Draw();
h02_entryZ_e->Draw("SAME");
h02_entryZ_n->Draw("SAME");
h02_entryZ_h->Draw("SAME");
h02_entryZ_mu->Draw("SAME");
c1->SetLogy();
```

<img src="images/image-20221104113829086.png" alt="image-20221104113829086" style="zoom:40%;" />

##### Photons and electrons momenta

```
h03_p_ph->Draw();
h03_p_e->Draw("SAME");
```

<img src="images/image-20221104113936955.png" alt="image-20221104113936955" style="zoom:40%;" />

##### Neutrons and other hadrons momenta

```
h04_p_n->Draw();
h04_p_h->Draw("SAME");
```

<img src="images/image-20221104114009703.png" alt="image-20221104114009703" style="zoom:40%;" />

##### Muons momenta:

```
h05_p_mu->Draw();
```

<img src="images/image-20221104114039529.png" alt="image-20221104114039529" style="zoom:40%;" />

##### Time of arrival at the detector entry point w.r.t. the bunch crossing

```
h06_T_ph->Draw();
h06_T_e->Draw("SAME");
h06_T_n->Draw("SAME");
h06_T_h->Draw("SAME");
h06_T_mu->Draw("SAME");
```

<img src="images/image-20221104114131508.png" alt="image-20221104114131508" style="zoom:40%;" />

#### The reconstructed detector hits of the BIB particles

##### Fill the histograms

```
root -l make_recoHitPlots_BIB.C+
```

##### Distributions of the tracker, calorimeters and muon detectors hits in the rho-z plane

```
h01_layout_trk->Draw("COLZ");
h01_layout_calo->Draw("COLZ,SAME");
h01_layout_muon->Draw("COLZ,SAME");
```

<img src="images/image-20221104114637906.png" alt="image-20221104114637906" style="zoom:40%;" />

##### Distributions of the vertex detector hits in the rho-z plane

```
h02_layout_VXD->Draw("COLZ");
```

<img src="images/image-20221104114750779.png" alt="image-20221104114750779" style="zoom:40%;" />

#####  Hit times w.r.t. the bunch crossing at the first layers of the vertex detector and the inner and outer trackers

```
h03_time_itB_L0->Draw();
h03_time_otB_L0->Draw("SAME");
h03_time_vxdB_L0->Draw("SAME");
```

<img src="images/image-20221104114913045.png" alt="image-20221104114913045" style="zoom:40%;" />

##### Hit energies deposited on the first layers of the vertex detector and the inner and outer trackers

```
h04_energy_vxdB_L0->Draw();
h04_energy_itB_L0->Draw("SAME");
h04_energy_otB_L0->Draw("SAME");
```

<img src="images/image-20221104115029857.png" alt="image-20221104115029857" style="zoom:40%;" />

##### ECAL barrel and ECAL+ endcap hit occupancies in the phi-theta and y-x planes

```
h05_occupancy_ecalB->Draw("COLZ");
h06_occupancy_ecalE->Draw("COLZ");
```

<img src="images/image-20221104115116993.png" alt="image-20221104115116993" style="zoom:40%;" />

<img src="images/image-20221104115134272.png" alt="image-20221104115134272" style="zoom:40%;" />

##### ECAL barrel and ECAL+ endcap hit occupancies in the rho-theta and y-z planes

```
h07_depth_ecalB->Draw("COLZ");
h08_depth_ecalE->Draw("COLZ");
```

<img src="images/image-20221104115219732.png" alt="image-20221104115219732" style="zoom:40%;" />

<img src="images/image-20221104115253198.png" alt="image-20221104115253198" style="zoom:40%;" />

##### HCAL barrel and HCAL+ endcap hit occupancies in the phi-theta and y-x planes

```
h09_occupancy_hcalB->Draw("COLZ");
h10_occupancy_hcalE->Draw("COLZ");
```

<img src="images/image-20221104115318940.png" alt="image-20221104115318940" style="zoom:40%;" />

<img src="images/image-20221104115334387.png" alt="image-20221104115334387" style="zoom:40%;" />

##### HCAL barrel and HCAL+ endcap hit occupancies in the rho-theta and y-z planes

```
h11_depth_hcalB->Draw("COLZ");
h12_depth_hcalE->Draw("COLZ");
```

<img src="images/image-20221104115428789.png" alt="image-20221104115428789" style="zoom:40%;" />

<img src="images/image-20221104115446478.png" alt="image-20221104115446478" style="zoom:40%;" />

##### Hit times at the ECAL and HCAL barrels

```
h13_time_ecalB->Draw();
h14_time_hcalB->Draw();
```

<img src="images/image-20221104115609671.png" style="zoom:40%;" />

<img src="images/image-20221104115630522.png" style="zoom:40%;" />

##### Hit energies depositid in the ECAL and HCAL barrels

```
h15_energy_ecalB->Draw();
c1->SetLogy();
h16_energy_hcalB->Draw();
```

<img src="images/image-20221104115727854.png" alt="image-20221104115727854" style="zoom:40%;" />

<img src="images/image-20221104115752461.png" alt="image-20221104115752461" style="zoom:40%;" />

## Back 1: Run the simulation remotely on `login.snowmass21.io`

After creating account by following the processes in Section 5. One can follow the tutorial from: https://confluence.infn.it/display/muoncollider/Tutorial

But if you want to modify the simulation file, I don't think run the simulation remotely is convenient.

## Backup 2: Some impractical way to install on macOS

### B2.1 Install manually

One can no longer install the Muon Collider Software by follow the instruction of https://confluence.infn.it/display/muoncollider/Installation. Because the DESY stops their SVN service in 2021, and this instruction depend on the SVN service of `svnsrv.desy.de` to download dependencies including `pathfinder` and `BBQ`.

### B2.2 Follow the docker setup tutorial on the official website of MuC Software

The official tutorial can not work on macOS. 
