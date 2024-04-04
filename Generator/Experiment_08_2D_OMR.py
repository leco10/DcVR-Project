"""
Leonard Constien 16.10.2023

social preference 2D VR experiment

restart kernel, remove variables!
"""
# imports:
from direct.showbase.ShowBase import ShowBase

from panda3d.core import loadPrcFile
loadPrcFile("VR_config.prc")

from panda3d.core import *

from direct.interval.LerpInterval import LerpColorInterval

# import task, a task is a procedure that is called every frame
from direct.task import Task

# install packages
import time
import numpy as np
import simplepbr
import random
import gltf
import os

#♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥
class VR(ShowBase):
    def __init__(self, cube_size =1024):
        super().__init__()
        
        # input
        self.fish_id = input('Enter fish ID: ')
        self.fish_age = input('Enter fish age (in months): ')
        self.fish_sex = input('Enter fish sex (m/f): ')               
        self.fish_length = input('Enter fish body length (in cm): ')
        
        # wait for uer input to start acquisition and experiment
        input('Press Enter to start the experiment, make sure tail tracking is running!') 
        
        # parameters
        self.habituation_duration = 120#600
        self.trial_duration = 40 # trial + intertrial
        self.trial_per_angle = 6
        self.alpha_list = []
        for i in range(self.trial_per_angle):
            alpha_group = [0,30,60,-30,-60]
            random.shuffle(alpha_group)
            self.alpha_list = self.alpha_list + alpha_group
        self.trial_number = self.trial_per_angle*5 # better even, exp time = hab + dur*(num)= 1800 s (with first habituation trial)
        self.max_inactive_trials = int((50*self.trial_number)/100)
        self.start = time.time()
        color = (0.05, 0.23, 0.85)
        self.v_max = 8 # cm/s231216_1
        self.av_max = 600
        self.conversion_factor = 0.8
        self.drift_velocity = -1.5*self.conversion_factor
        self.col_num = self.trial_duration*60 + 100 # 60 fps, +100 for saftey
        self.mold_height = 1.5#2.25
        self.projector_distance = -29
        self.projector_angle =98.65#98.45
        self.offset = 4.62
        self.cam_fov = 20.5#18.631
        self.starting_height = 0.6
        
        # variables
        self.end_of_trial = False
        self.experiment_running = False
        self.drift_0 = False
        self.drift_non0 = False
        self.drift_vector = Vec3(0,1,0)
        self.timer = 0
        self.counter = 0
        self.inactive_trial_counter = 0
        self.velocity = 0
        self.actual_drift_velocity = 0
        self.T = [0,0]
        self.trial_vel = []
        self.trial_ang_vel = []
        self.trial_x_positions = [0,0]
        self.trial_y_positions = [0,0]
        self.trial_headings = [0,0]
        self.trial_starts = []
        self.inter_trial_starts = []
        self.X_POS=[[0]*self.col_num]
        self.Y_POS=[[0]*self.col_num]
        self.HEAD=[[0]*self.col_num]
        self.timestamp = []
        
        # create new window, render and camera and create polygon for fading
        # window properties for new window
        wp = WindowProperties()
        wp.setSize(720,720)
        wp.setOrigin(1920,0) #1920,0
        
        # open new window for observing the bowl, set the gsg to main window to transfer textures
        self.bowl_win = self.openWindow(props=wp,gsg=self.win)
        # create render and camera
        self.bowl_render = NodePath('render2')  # the string parameter is important
        self.bowl_cam = self.camList[1] # is the camera of self.bowl_win
        # make polygon with Cardmaker for fading
        self.card = CardMaker("fade card")
        self.card.setFrame(-200,200,-200,200)
        self.cardNP = NodePath(self.card.generate())
        self.cardNP.reparentTo(self.bowl_cam)
        self.cardNP.setPos(0,5,0)
        self.cardNP.setTransparency(TransparencyAttrib.MAlpha)
        self.cardNP.setColor((0,0,0,0))
        
        # create environment in main window
        self.floor = gltf._loader.load_model(r"E:\leos data\VR environments\pebbles.glb")
        self.floor = NodePath(self.floor) #convert modelroot to nodepath
        self.floor.reparentTo(self.render)
        
        # ambient lightning
        alight2 = AmbientLight('alight')
        alight2.setColor((0.7,0.7,0.7,1)) #for naturalistic environment
        
        alnp2 = self.render.attachNewNode(alight2)
        self.render.setLight(alnp2)
        
        # set fog
        expfog = Fog("Scene-wide exponential Fog object")
        expfog.setColor(*color)
        expfog.setExpDensity(0.03)
        self.floor.setFog(expfog)
        
        # create subject, attach avatar and adjust overview camera
        self.subject = NodePath('subject')
        self.subject.reparentTo(self.render)
        self.subject.setPos(0,0,self.starting_height)
        # create avatar
        self.avatar = self.loader.load_model("/e/leos data/VR environments/Dc_animation_low_poly.egg")
        self.avatar = NodePath(self.avatar)
        self.avatar.setScale(2.2)
        self.avatar.setColor((1,0,0,0))
        self.avatar.setHpr(180,0,0)
        self.avatar.reparentTo(self.subject)
        
        # adjust overview camera
        self.cam.reparentTo(self.subject)
        self.cam.setPos(0,-5,8)
        self.cam.lookAt(self.avatar)
        self.camLens.setFov(75)
        self.camLens.setNearFar(1/100,100)
        
        # create bowl in bowl window, adjust bowl camera and setup projector from subject position to bowl
        self.bowl = gltf._loader.load_model(r"E:\leos data\VR environments\flat150_bowl.glb")
        self.bowl = NodePath(self.bowl) #convert modelroot to nodepath
        self.bowl.reparentTo(self.bowl_render)
        # adjust camera
        self.bowl_cam.reparentTo(self.bowl)
        self.bowl_cam.node().getLens().setNearFar(1/100,100)
        self.bowl_cam.setPos(0,self.offset,self.projector_distance) 
        self.bowl_cam.setHpr(0,self.projector_angle,0) 
        self.bowl_cam.node().getLens().setFov(self.cam_fov)
        # create projector to project cube map onto bowl    
        self.bowl_proj = self.bowl.attachNewNode(LensNode('proj'))
        self.lens = PerspectiveLens()
        self.lens.setNearFar(1/1000, 10000)
        self.lens.setFov(60) #needs to be 60!
        self.bowl_proj.node().setLens(self.lens)
        #self.bowl_proj.node().showFrustum()
        self.bowl_proj.setPos(0,0, self.mold_height)
        self.bowl_proj.setHpr(0,90,0)
        
        # setup cube cameras in main window, generate cube map and back-project to bowl in bowl window
        # create cube map cameras, and move them with the subject
        self.cube_cams = self.render.attachNewNode('cube_cams')
        # create camera mask to hide avatar from cube cams
        self.CubeCameraMask = BitMask32.bit(1)
        self.avatar.hide(self.CubeCameraMask)
        # create camera mask for cube cams
        cube_map_buffer = self.win.makeCubeMap('cube_map', cube_size, self.cube_cams,self.CubeCameraMask)
        cube_map_buffer.set_inverted(True)
        [icam.node().getLens().setNearFar(1/1000,100) for icam in self.cube_cams.getChildren()]
        # back-project cube map onto bowl in other window from the location of the subject
        self.bowl.setTexGen(TextureStage.getDefault(), TexGenAttrib.MWorldPosition)
        self.bowl.projectTexture(TextureStage.getDefault(),cube_map_buffer.getTexture(), self.bowl_proj)
        self.bowl.setTexScale(TextureStage.getDefault(), 1)
        self.bowl.setTexHpr(TextureStage.getDefault(), 0,0,0)
        self.bowl.setTexPos(TextureStage.getDefault(), -0.50,-0.5,0)
        
        # simplepbr for main window
        simplepbr.init(enable_fog=True, enable_shadows=True,render_node=self.floor)
        # for bowl window
        simplepbr.init(enable_fog=True,enable_shadows=True,window=self.bowl_win,camera_node=self.bowl_cam)
        
        self.exp_start = 'Experiment commenced at ' + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print(self.exp_start)
        print('Drift angle list: ' + str(self.alpha_list))
        print("Habituation running for " + str(self.habituation_duration/60) + " min")
        self.end = False
        self.start_of_experiment = False
        
        # tasks
        self.taskMgr.add(self.protocol, "protocol")
        self.taskMgr.add(self.move_sub, "move_sub")
        
    def protocol(self, task):
        if self.counter <= self.trial_number:
            #initial fade in
            if self.timer == 0 and self.counter == 0:
                LerpColorInterval(self.cardNP, 7, (0,0,0,0), startColor=(0,0,0,1)).start()
               
            # exp time
            self.t = time.time()-self.start
            # get delta t
            self.T.append(self.t)
            self.dt = self.T[-1]-self.T[-2]
            # progress timer
            self.timer += self.dt
            
            #end of habituation
            if self.timer > self.habituation_duration or self.experiment_running == True:
                
                # save trial data
                self.trial_x_positions.append(self.subject.getPos().x)
                self.trial_y_positions.append(self.subject.getPos().y)
                self.trial_headings.append(self.subject.getHpr().x)
                
                #start of trial 
                if self.timer > self.trial_duration or self.counter == 0:
                    #experiment is commencing
                    if self.counter == 0:
                        self.experiment_running = True
                         
                    # check if it was an inactive trial and save data
                    if self.counter > 0:
                       #self.check_exclusion()
                   
                       # fill up position and heading arrays with 0s
                       self.trial_x_positions = self.trial_x_positions + [0]*(self.col_num-len(self.trial_x_positions))
                       self.trial_y_positions = self.trial_y_positions + [0]*(self.col_num-len(self.trial_y_positions))
                       self.trial_headings = self.trial_headings + [0]*(self.col_num-len(self.trial_headings))
                   
                       # save positions in matrix
                       self.X_POS = np.append(self.X_POS,[self.trial_x_positions],axis=0)
                       self.Y_POS = np.append(self.Y_POS,[self.trial_y_positions],axis=0)
                       self.HEAD = np.append(self.HEAD,[self.trial_headings],axis=0)
                      
                       # reset trial variables
                       self.trial_x_positions = [0,0]
                       self.trial_y_positions = [0,0]
                       self.trial_headings = [0,0]
                       self.trial_vel = []
                       self.trial_ang_vel = []
                   
                    # set alpha, drift vector and start lerping floor hpr
                    if self.counter < self.trial_number:
                        self.alpha = self.alpha_list[self.counter]
                        print("Starting trial " + str(self.counter+1) + "/" + str(self.trial_number) + ", drift angle: " + str(self.alpha))
                    quat = self.subject.getQuat()
                    forwardVec = quat.getForward()
                    self.drift_vector = Vec3(forwardVec.x*np.cos(self.alpha*np.pi/180) - forwardVec.y*np.sin(self.alpha*np.pi/180), forwardVec.x*np.sin(self.alpha*np.pi/180) +forwardVec.y*np.cos(self.alpha*np.pi/180),0)
                    
                    # save trial starts
                    self.trial_starts.append(self.timestamp[-1])
                    
                    # reset other stuff and progress counter
                    self.timer = 0
                    self.drift_0 = False
                    self.drift_non0 = False
                    self.counter += 1
                    
                elif self.timer > self.trial_duration/2 and self.drift_0==False:
                    print("Inter-trial")
                    # save inter-trial starts
                    self.inter_trial_starts.append(self.timestamp[-1])
                    self.drift_0 = True
                    self.actual_drift_velocity = 0
                    
                # activate it else    
                elif self.drift_non0 == False:
                    self.drift_non0 = True
                    self.actual_drift_velocity = self.drift_velocity
                        
            return Task.cont   
        
        elif self.end == False:
            time.sleep(1)
            print('Experiment finished at ' + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
            self.end = True
            self.end_experiment()
        #---------------------------------------------------------------------- 
    def move_sub (self, task):
        
        # get vigor and convert to velocity if above threshold
        print(end='\x1b[2K')
        with open(r'C:\Users\jlab\Desktop\Leo\TailTracker\vigor_data.txt','r') as f:
            str_input = f.readline()
            result = [val for val in str_input.split(',')]
   
        if not '' in result:
            vdata = float(result[0]) #cm/s
            avdata = float(result[1]) #degree/s
            if vdata <= self.v_max: # check if below maximal velocity
                self.velocity = vdata*self.conversion_factor
            else:
                self.velocity = self.v_max*self.conversion_factor
            if avdata <= self.av_max:
                self.angular_velocity = avdata
            else:
                self.angular_velocity = self.av_max
            
            self.trial_vel.append(self.velocity)
            self.trial_ang_vel.append(self.angular_velocity)
            self.timestamp.append(float(result[2]))
        
        #get subject Hpr
        sub_hpr = self.subject.getHpr()
        sub_pos = self.subject.getPos()
        
        quat = self.subject.getQuat()
        forwardVec = quat.getForward()
        
        #change position
        final_velocity = self.velocity*self.dt 
        final_drift = self.actual_drift_velocity*self.dt
        
        self.floor.setPos(self.floor.getPos() - self.drift_vector*final_drift)
        
        self.subject.setPos(self.subject.getPos() + forwardVec*final_velocity)
        self.cube_cams.setPos(self.subject.getPos())# update cube cam position, not sure why attaching the node doesn't work

        final_angular_velocity = self.angular_velocity*self.dt
        sub_hpr.x += final_angular_velocity
        self.subject.setHpr(sub_hpr)         
        self.cube_cams.setHpr(sub_hpr) # update cube cam hpr, not sure why attaching the node doesn't work
        
        #teleport
        if sub_pos.y < self.floor.getPos().y-240:
            self.subject.setPos(Vec3(self.floor.getPos().x,self.floor.getPos().y,self.starting_height) + forwardVec*final_velocity)
           
        if sub_pos.y >self.floor.getPos().y+240:
            self.subject.setPos(Vec3(self.floor.getPos().x,self.floor.getPos().y,self.starting_height) + forwardVec*final_velocity)
           
        if sub_pos.x < self.floor.getPos().x-240:
            self.subject.setPos(Vec3(self.floor.getPos().x,self.floor.getPos().y,self.starting_height) + forwardVec*final_velocity)
        
        if sub_pos.x > self.floor.getPos().x+240:
            self.subject.setPos(Vec3(self.floor.getPos().x,self.floor.getPos().y,self.starting_height) + forwardVec*final_velocity)
        
        return Task.cont
    
    def check_exclusion(self):
        if any(self.trial_vel) == False:
            if any(self.trial_ang_vel) == False:
                self.inactive_trial_counter += 1
                print(str(self.inactive_trial_counter) + "/" + str(self.max_inactive_trials) + " inactive trials!")
                
        if self.inactive_trial_counter == self.max_inactive_trials:
            print("Experiment stopped! Fish is inactive!")
            # save to file
            #np.savetxt('E:\leos data\Exp07_social_preference\Inactive_' + str(self.fish_id) + '.txt',[0])
            VR_instance.destroy()
            time.sleep(1000)
        
    def end_experiment(self):
        
        # save to file
        path = os.path.join('E:\leos data\Exp08_2D_OMR',str(self.fish_id)).replace("\\","/")
        os.mkdir(path) 
        np.savetxt(os.path.join(path,'X_POS_' + str(self.fish_id) + '.txt').replace("\\","/"), self.X_POS,'%s')
        np.savetxt(os.path.join(path,'Y_POS_' + str(self.fish_id) + '.txt').replace("\\","/"), self.Y_POS,'%s')
        np.savetxt(os.path.join(path,'HEAD_' + str(self.fish_id) + '.txt').replace("\\","/"),  self.HEAD,'%s')
        np.savetxt(os.path.join(path,'MISC' + str(self.fish_id) + '.txt').replace("\\","/"),  self.T,'%s',header='Meta data: \n' + 'Fish ID: ' + str(self.fish_id) + '\n' + 'Fish age (months): ' + self.fish_age + '\n' + 'Fish sex: '
                   + self.fish_sex + '\n' + 'Fish length (cm): ' + self.fish_length + "\n" + 'alpha list: ' + str(self.alpha_list) + "\n"  + 'trials per angle: ' + str(self.trial_per_angle) + "\n" + 'trial starts: ' + str(self.trial_starts) + "\n" + "inter-trial starts: " + str(self.inter_trial_starts) +
                   '\n' + 'Time: ')

        # terminate program and clean up
        VR_instance.destroy()

     #♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥
    
VR_instance = VR() #read out variables with VR_instance.variable
VR_instance.run()  


