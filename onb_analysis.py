import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.transform import Rotation as R, RotationSpline
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os

def fname_to_dt(df:pd.DataFrame)->pd.DataFrame:
    df_meta = pd.DataFrame([{"DateTime":s[0]+s[1],
    "plant label":s[2],
    "lamina no.":int(s[3][-3:])} for s in df.file.str.replace(".txt","").str.split("_")])
    df_meta["DateTime"] = pd.to_datetime(df_meta["DateTime"])
    df_meta_timedelta = df_meta["DateTime"].diff()
    df_meta_timedelta.iloc[0] = pd.to_timedelta(0)
    df_meta["TimeDelta"] = df_meta_timedelta
    df_meta["TimeDelta(min)"] = df_meta["TimeDelta"]/pd.Timedelta(minutes=1)
    df_meta["Duration"] = df_meta_timedelta.cumsum()
    df_meta["Duration(min)"] = df_meta["Duration"]/pd.Timedelta(minutes=1)
    df = pd.concat([df.reset_index(drop=True),df_meta],axis=1)
    return df

def p_to_xz_plane(df:pd.DataFrame)->pd.DataFrame:
    p_azim = np.arctan2(df["p_y"],df["p_x"]).values

    labels = ["p","basal_vec","pca1","pca2","pca3"]
    labels_xyz = []
    for a in labels:
        labels_xyz += [f"{a}_x2",f"{a}_y2",f"{a}_z2"]

    rotated_vecs = []

    for i in range(len(p_azim)):
        rot = R.from_rotvec(np.array([0,0,1])*(-p_azim[i]))
        p2 = rot.apply(df[["p_x","p_y","p_z"]].values[i])
        for a in labels[1:]:
            p2 = np.hstack([p2,rot.apply(df[[f"{a}_x",f"{a}_y",f"{a}_z"]].values[i])])
        rotated_vecs.append(p2)
    df_rotated_vecs = pd.DataFrame(rotated_vecs,columns=labels_xyz)
    return df_rotated_vecs

def draw_veccoord_plot(df:pd.DataFrame,fname=""):
    fig,axes = plt.subplots(figsize=(24,4),ncols=5,nrows=1)
    plt.subplots_adjust(wspace=0.3)
    for n,variable in enumerate(["p","basal_vec","pca1","pca2","pca3"]):
        ax = axes[n]
        for a in ["x2","y2","z2"]:
            for i,f in enumerate([sns.lineplot,sns.scatterplot]):
                if i==0:
                    l=f"{variable}_{a}"
                else:
                    l=None
                f(
                    data=df,
                    x="Duration(min)",
                    y=f"{variable}_{a}",
                    label=l,
                    ax=ax
                )
        ax.legend()
        ax.set_ylim(-1.1,1.1)
        if len(fname)>0:
            fig.savefig(f"graph/{fname}.png",dpi=250,bbox_inches="tight")
    plt.show()

def calc_add_bvw(df:pd.DataFrame)->pd.DataFrame:
    variable = "basal_vec"
    df_1 = df[[f"{variable}_{a}" for a in ["x2","y2","z2"]]]
    variable = "pca2"
    df_2 = df[[f"{variable}_{a}" for a in ["x2","y2","z2"]]]
    w_list = []
    for i in range(len(df_1)):
        w_list.append(np.cross(df_1.iloc[i,:].values,df_2.iloc[i,:].values))
    df["bv_w_x2"] = np.array(w_list)[:,0]
    df["bv_w_y2"] = np.array(w_list)[:,1]
    df["bv_w_z2"] = np.array(w_list)[:,2]
    return df

class LeafMovRot():
    def __init__(self):
        self.onb_cols = ["basal_vec_x2","basal_vec_y2","basal_vec_z2",
                "pca2_x2","pca2_y2","pca2_z2",
                "bv_w_x2","bv_w_y2","bv_w_z2"]
        
    def fit(self,df:pd.DataFrame):
        self.df = df
        onb_list = []
        for i in range(len(df)):
            onb = df.loc[i,self.onb_cols].values.reshape(3,3).T
            onb_list.append(onb)
        self.rotations = R.from_matrix(onb_list) 
        self.t_values = df["Duration(min)"]
        spline = RotationSpline(self.t_values, self.rotations)

        self.t_interp = np.arange(0, self.t_values.max(), 1)
        self.R_interp = spline(self.t_interp)
        self.interp_mtx = self.R_interp.as_matrix()
        self.interp_rotvec = self.R_interp.as_rotvec()

        rot_relative = self.rotations[0].inv()*self.rotations
        self.rot_relative_pm = rot_relative.as_rotvec()

        interp_relative = self.R_interp[0].inv()*self.R_interp
        self.interp_relative_pm = interp_relative.as_rotvec()

        t_diffs= self.t_values.diff().values[1:]
        rot_diffs = self.rotations[:-1].inv()*self.rotations[1:]
        rot_diffs_pm = (rot_diffs.as_rotvec().T/t_diffs).T
        self.rot_diffs_pm_lamaxis = self.rotations[:-1].inv().apply(rot_diffs_pm)

        self.delta_R = self.R_interp[:-1].inv() * self.R_interp[1:]
        self.delta_R_rev =  self.R_interp[1:] * self.R_interp[:-1].inv()
        self.omega_world = self.delta_R.as_rotvec()
        #self.omega_world_rev = self.delta_R_rev.as_rotvec()
        self.omega_body = self.R_interp[:-1].inv().apply(self.omega_world)

class LeafMovRot2():
    def __init__(self):
        self.onb_cols = ["pca1_x2","pca1_y2","pca1_z2",
                "pca2_x2","pca2_y2","pca2_z2",
                "pca3_x2","pca3_y2","pca3_z2"]
        
    def fit(self,df:pd.DataFrame):
        self.df = df
        onb_list = []
        for i in range(len(df)):
            onb = df.loc[i,self.onb_cols].values.reshape(3,3).T
            onb_list.append(onb)
        self.rotations = R.from_matrix(onb_list) 
        self.t_values = df["Duration(min)"]
        spline = RotationSpline(self.t_values, self.rotations)

        self.t_interp = np.arange(0, self.t_values.max(), 1)
        self.R_interp = spline(self.t_interp)
        self.interp_mtx = self.R_interp.as_matrix()
        self.interp_rotvec = self.R_interp.as_rotvec()

        rot_relative = self.rotations[0].inv()*self.rotations
        self.rot_relative_pm = rot_relative.as_rotvec()

        interp_relative = self.R_interp[0].inv()*self.R_interp
        self.interp_relative_pm = interp_relative.as_rotvec()

        t_diffs= self.t_values.diff().values[1:]
        rot_diffs = self.rotations[:-1].inv()*self.rotations[1:]
        rot_diffs_pm = (rot_diffs.as_rotvec().T/t_diffs).T
        self.rot_diffs_pm_lamaxis = self.rotations[:-1].inv().apply(rot_diffs_pm)

        #modified at 20250819
        # 相対回転（body基準）: A1^{-1} A2
        self.delta_R    = self.R_interp[:-1].inv() * self.R_interp[1:]
        # 相対回転（world基準）: A2 A1^{-1}
        self.delta_R_w  = self.R_interp[1:] * self.R_interp[:-1].inv()

        # 回転ベクトル（= 角速度 × Δt）。Δt=1 なので、この値がそのまま「1ステップ当たりの角速度成分」に等しい。
        self.omega_body  = self.delta_R.as_rotvec()          # body基準
        self.omega_world = self.R_interp[:-1].apply(self.omega_body)
        # あるいは world を直接：
        # self.omega_world = self.delta_R_w.as_rotvec()


def draw_series_onb(lmr:LeafMovRot,savedir="",elev=30, azim=45):
    if len(savedir)>0:
        savedir2 = "graph/"+savedir
        os.makedirs(savedir2,exist_ok=True)
    
    colors = ["blue", "orange", "indigo"]
    labels = ["basal_vec", "pca2", "bv_w"]
    labels2 = ["PD(basal)", "ML", "Ad(basal)"]

    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 14

    px,py,pz = lmr.df[["p_x2","p_y2","p_z2"]].mean(axis=0).values
    omega_world = np.vstack([lmr.omega_world,lmr.omega_world[-1,:]])
    #omega_world_rev = np.vstack([lmr.omega_world_rev,lmr.omega_world_rev[-1,:]])

    #for n in range(1):
    for n in range(lmr.interp_mtx.shape[0]):
        ow = omega_world[n,:]*100
        #ow_rev = omega_world_rev[n,:]*100
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # 球面の描画（透明）
        u, v = np.meshgrid(np.linspace(0, 2*np.pi, 50),
                        np.linspace(0, np.pi, 50))
        sphere_x = np.cos(u) * np.sin(v)
        sphere_y = np.sin(u) * np.sin(v)
        sphere_z = np.cos(v)
        ax.plot_surface(sphere_x, sphere_y, sphere_z, color='lightblue', alpha=0.1, edgecolor='none')

        # 原点と2ベクトルで三角形を構成
        v0 = np.array([0, 0, 0])
        v1 = np.array([lmr.interp_mtx[n, 0, 0], lmr.interp_mtx[n, 1, 0], lmr.interp_mtx[n, 2, 0]])
        v2 = np.array([(lmr.interp_mtx[n, 0, 0]+lmr.interp_mtx[n, 0, 1])/2,
                    (lmr.interp_mtx[n, 1, 0]+lmr.interp_mtx[n, 1, 1])/2,
                    (lmr.interp_mtx[n, 2, 0]+lmr.interp_mtx[n, 2, 1])/2])
        v3 = np.array([(lmr.interp_mtx[n, 0, 0]-lmr.interp_mtx[n, 0, 1])/2,
                    (lmr.interp_mtx[n, 1, 0]-lmr.interp_mtx[n, 1, 1])/2,
                    (lmr.interp_mtx[n, 2, 0]-lmr.interp_mtx[n, 2, 1])/2])

        triangle = [v0, v3, v1, v2]
        triangle_poly = Poly3DCollection([triangle], alpha=0.5, color='darkgreen')
        ax.add_collection3d(triangle_poly)

        ax.plot(
                    [lmr.interp_mtx[n, 0, 0]/2,(lmr.interp_mtx[n, 0, 0]+lmr.interp_mtx[n, 0, 1])/2],
                    [lmr.interp_mtx[n, 1, 0]/2,(lmr.interp_mtx[n, 1, 0]+lmr.interp_mtx[n, 1, 1])/2],
                    [lmr.interp_mtx[n, 2, 0]/2,(lmr.interp_mtx[n, 2, 0]+lmr.interp_mtx[n, 2, 1])/2],
                    color='darkgreen',
                    linewidth=2
                )

        for i2 in range(3):
            ax.plot(
                        [0,lmr.interp_mtx[n, 0, i2]],
                        [0,lmr.interp_mtx[n, 1, i2]],
                        [0,lmr.interp_mtx[n, 2, i2]],
                        color='black',
                        linestyle='dotted',
                        linewidth=2
                    )
        
        ax.plot(
                    [-px,0],
                    [-py,0],
                    [-pz,0],
                    color='darkolivegreen',
                    linewidth=4,
                    alpha=0.6
                )
        
        ax.plot(
                    [0,ow[0]],
                    [0,ow[1]],
                    [0,ow[2]],
                    color='salmon',
                    linewidth=2,
                    alpha=0.8
                )

        # 軌跡・点・線の描画ループ
        for i, l in enumerate(labels):
            x = lmr.df[f"{l}_x2"]
            y = lmr.df[f"{l}_y2"]
            z = lmr.df[f"{l}_z2"]
            x_interp = lmr.interp_mtx[:, 0, i]
            y_interp = lmr.interp_mtx[:, 1, i]
            z_interp = lmr.interp_mtx[:, 2, i]

            # 元データ点
            ax.scatter(x, y, z, color='black', s=20)

            # 補間軌跡
            ax.plot(x_interp, y_interp, z_interp, color=colors[i], label=labels2[i], linewidth=2)

        # 軸・ラベル
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"Duration Time: {n:04d} min") #全てOK
        #ax.set_title(f"Duration Time: {n:04d} min: {np.dot(interp_mtx[n, :, 0],interp_mtx[n, :, 1]):.03f} {np.dot(interp_mtx[n, :, 0],interp_mtx[n, :, 2]):.03f}")

        # 軸のスケールを等しく
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])
        ax.set_zlim([-1.1, 1.1])
        ax.set_box_aspect([1, 1, 1])  # matplotlib>=3.3

        # 視点角度の指定（elevation, azimuth）
        ax.view_init(elev=elev, azim=azim)

        # 凡例
        ax.legend()

        if len(savedir)>0:
            fig.savefig(f"{savedir2}/{n:05d}.png",dpi=150,bbox_inches="tight")
        plt.close()

def calc_grav_vec(df):
        body_grav_vec_list = []

        for n in range(len(df)):
                rot_mtx = np.array([df.iloc[n,:].filter(like="x2").values[2:5],
                        df.iloc[n,:].filter(like="y2").values[2:5],
                        df.iloc[n,:].filter(like="z2").values[2:5]])

                rot = R.from_matrix(rot_mtx)
                #rot.inv().apply(rot_mtx.T[0]).astype(float).round(10)
                body_grav_vec = rot.inv().apply(np.array([0,0,-1],dtype=np.float64))
                body_grav_vec_list.append(body_grav_vec)

        df_grav = pd.DataFrame(body_grav_vec_list,
                        columns=["grav_vec_x_body",
                                                "grav_vec_y_body",
                                                "grav_vec_z_body"])
        df_grav["Duration(min)"] = df["Duration(min)"].values

        return df_grav