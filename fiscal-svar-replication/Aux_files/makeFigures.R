#figure theme
fig_thm = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background=element_blank(),
                panel.border =element_rect(colour='black',fill=NA),
                axis.text=element_text(color='black'),
                axis.ticks=element_line(color='black'),
                plot.title = element_text(face="bold",hjust=0.5),
                text=element_text(family = "Times"),
                legend.title=element_blank(),
                legend.key = element_blank())

#plot shocks (Figure A1)
eps.df = data.frame(G = eps[,"G"], T = eps[,"T"], Y= eps[,"Y"],
                    G.iv = eps.iv[,"G"], T.iv = eps.iv[,"T"], Y.iv = eps.iv[,"Y"])

#government spending shocks (Figure A1, top panel)
fig.A1.1 = ggplot(data=eps.df,aes(x=as.Date(rownames(eps)),group=1))+
  geom_line(aes(y=G.iv,color="'black'"))+
  geom_line(aes(y=G,color="'red'"))+
  scale_x_date(expand=c(0,0))+
  labs(x="Year", y="shocks", title="Government spending")+
  scale_color_manual(values = c("black","red"))+
  theme(legend.position = "none")+
  fig_thm

#tax shocks (Figure A1, middle panel)
fig.A1.2 = ggplot(data=eps.df,aes(x=as.Date(rownames(eps)),group=1))+
  geom_line(aes(y=T.iv,color="'black'"))+
  geom_line(aes(y=T,color="'red'"))+
  scale_x_date(expand=c(0,0))+
  labs(x="Year", y="shocks", title="Tax Revenue")+
  scale_color_manual(values = c("black","red"))+
  theme(legend.position = "none")+
  fig_thm

#output shocks (Figure A1, bottom panel)
fig.A1.3 = ggplot(data=eps.df,aes(x=as.Date(rownames(eps)),group=1))+
  geom_line(aes(y=Y.iv,color="'black'"))+
  geom_line(aes(y=Y,color="'red'"))+
  scale_x_date(expand=c(0,0))+
  labs(x="Year", y="shocks", title="Output")+
  scale_color_manual(values = c("black","red"))+
  theme(legend.position = "none")+
  fig_thm

fig.A1 = plot_grid(fig.A1.1,fig.A1.2,fig.A1.3,ncol=1)

#plot multipliers (Figure 1)

multipliers.df = data.frame(horizon = seq(0,size(MULT,1)-1),
                            m.G = MULT[,"G"],
                            m.G.lower68 = MULT[,"G"]-se.mult[,"G"],
                            m.G.upper68 = MULT[,"G"]+se.mult[,"G"],
                            m.G.iv = MULT.iv[,"G"],
                            m.G.lower68.iv = MULT.iv[,"G"]-MULT.se.iv[,"G"],
                            m.G.upper68.iv = MULT.iv[,"G"]+MULT.se.iv[,"G"],
                            m.T = MULT[,"T"],
                            m.T.lower68 = MULT[,"T"]-se.mult[,"T"],
                            m.T.upper68 = MULT[,"T"]+se.mult[,"T"],
                            m.T.iv = MULT.iv[,"T"],
                            m.T.lower68.iv = MULT.iv[,"T"]-MULT.se.iv[,"T"],
                            m.T.upper68.iv = MULT.iv[,"T"]+MULT.se.iv[,"T"],
                            dm.G = DMULT[,"G"],
                            dm.G.lower68 = DMULT[,"G"]-se.dmult[,"G"],
                            dm.G.upper68 = DMULT[,"G"]+se.dmult[,"G"],
                            dm.G.iv = DMULT.iv[,"G"],
                            dm.G.lower68.iv = DMULT.iv[,"G"]-DMULT.se.iv[,"G"],
                            dm.G.upper68.iv = DMULT.iv[,"G"]+DMULT.se.iv[,"G"],
                            dm.T = DMULT[,"T"],
                            dm.T.lower68 = DMULT[,"T"]-se.dmult[,"T"],
                            dm.T.upper68 = DMULT[,"T"]+se.dmult[,"T"],
                            dm.T.iv = DMULT.iv[,"T"],
                            dm.T.lower68.iv = DMULT.iv[,"T"]-DMULT.se.iv[,"T"],
                            dm.T.upper68.iv = DMULT.iv[,"T"]+DMULT.se.iv[,"T"]
                            )

fig.1.1 = ggplot(data=multipliers.df,aes(x=horizon),group=1)+
  geom_line(aes(y=dm.G.iv),color="black")+
  geom_line(aes(y=dm.G),color="red")+
  geom_line(aes(y=dm.G.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=dm.G.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=dm.G.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=dm.G.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="Dynamic: g-shock")+
  theme(legend.position = "none")+
  fig_thm

fig.1.2 = ggplot(data=multipliers.df,aes(x=horizon),group=1)+
  geom_line(aes(y=dm.T.iv),color="black")+
  geom_line(aes(y=dm.T),color="red")+
  geom_line(aes(y=dm.T.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=dm.T.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=dm.T.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=dm.T.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="Dynamic: tax shock")+
  theme(legend.position = "none")+
  fig_thm

fig.1.3 = ggplot(data=multipliers.df,aes(x=horizon),group=1)+
  geom_line(aes(y=m.G.iv),color="black")+
  geom_line(aes(y=m.G),color="red")+
  geom_line(aes(y=m.G.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=m.G.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=m.G.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=m.G.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="Cumulative: g shock")+
  theme(legend.position = "none")+
  fig_thm

fig.1.4 = ggplot(data=multipliers.df,aes(x=horizon),group=1)+
  geom_line(aes(y=m.T.iv),color="black")+
  geom_line(aes(y=m.T),color="red")+
  geom_line(aes(y=m.T.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=m.T.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=m.T.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=m.T.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="Cumulative: tax shock")+
  theme(legend.position = "none")+
  fig_thm

fig.1 = plot_grid(fig.1.1,fig.1.2,fig.1.3,fig.1.4,ncol=2)

#Figure A2 (impulse response functions)

irfs.df = data.frame(horizon = seq(0,size(IRF,3)-1),
                     #government spending: over identified
                     G.to.G = IRF["G","G",],
                     G.to.T = IRF["G","T",],
                     G.to.Y = IRF["G","Y",],
                     G.to.G.lower68 = IRF["G","G",]-se.irf[,"GG"],
                     G.to.T.lower68 = IRF["G","T",]-se.irf[,"GT"],
                     G.to.Y.lower68 = IRF["G","Y",]-se.irf[,"GY"],
                     G.to.G.upper68 = IRF["G","G",]+se.irf[,"GG"],
                     G.to.T.upper68 = IRF["G","T",]+se.irf[,"GT"],
                     G.to.Y.upper68 = IRF["G","Y",]+se.irf[,"GY"],
                     #government spending: just identified
                     G.to.G.iv = IRF.iv["G","G",],
                     G.to.T.iv = IRF.iv["G","T",],
                     G.to.Y.iv = IRF.iv["G","Y",],
                     G.to.G.lower68.iv = IRF.iv["G","G",]-IRF.se.iv[,"GG"],
                     G.to.T.lower68.iv = IRF.iv["G","T",]-IRF.se.iv[,"GT"],
                     G.to.Y.lower68.iv = IRF.iv["G","Y",]-IRF.se.iv[,"GY"],
                     G.to.G.upper68.iv = IRF.iv["G","G",]+IRF.se.iv[,"GG"],
                     G.to.T.upper68.iv = IRF.iv["G","T",]+IRF.se.iv[,"GT"],
                     G.to.Y.upper68.iv = IRF.iv["G","Y",]+IRF.se.iv[,"GY"],
                     #tax revenue: over identified
                     T.to.G = IRF["T","G",],
                     T.to.T = IRF["T","T",],
                     T.to.Y = IRF["T","Y",],
                     T.to.G.lower68 = IRF["T","G",]-se.irf[,"TG"],
                     T.to.T.lower68 = IRF["T","T",]-se.irf[,"TT"],
                     T.to.Y.lower68 = IRF["T","Y",]-se.irf[,"TY"],
                     T.to.G.upper68 = IRF["T","G",]+se.irf[,"TG"],
                     T.to.T.upper68 = IRF["T","T",]+se.irf[,"TT"],
                     T.to.Y.upper68 = IRF["T","Y",]+se.irf[,"TY"],
                     #tax revenue: just identified
                     T.to.G.iv = IRF.iv["T","G",],
                     T.to.T.iv = IRF.iv["T","T",],
                     T.to.Y.iv = IRF.iv["T","Y",],
                     T.to.G.lower68.iv = IRF.iv["T","G",]-IRF.se.iv[,"TG"],
                     T.to.T.lower68.iv = IRF.iv["T","T",]-IRF.se.iv[,"TT"],
                     T.to.Y.lower68.iv = IRF.iv["T","Y",]-IRF.se.iv[,"TY"],
                     T.to.G.upper68.iv = IRF.iv["T","G",]+IRF.se.iv[,"TG"],
                     T.to.T.upper68.iv = IRF.iv["T","T",]+IRF.se.iv[,"TT"],
                     T.to.Y.upper68.iv = IRF.iv["T","Y",]+IRF.se.iv[,"TY"],
                     #output: over identified
                     Y.to.G = IRF["Y","G",],
                     Y.to.T = IRF["Y","T",],
                     Y.to.Y = IRF["Y","Y",],
                     Y.to.G.lower68 = IRF["Y","G",]-se.irf[,"YG"],
                     Y.to.T.lower68 = IRF["Y","T",]-se.irf[,"YT"],
                     Y.to.Y.lower68 = IRF["Y","Y",]-se.irf[,"YY"],
                     Y.to.G.upper68 = IRF["Y","G",]+se.irf[,"YG"],
                     Y.to.T.upper68 = IRF["Y","T",]+se.irf[,"YT"],
                     Y.to.Y.upper68 = IRF["Y","Y",]+se.irf[,"YY"],
                     #tax revenue: just identified
                     Y.to.G.iv = IRF.iv["Y","G",],
                     Y.to.T.iv = IRF.iv["Y","T",],
                     Y.to.Y.iv = IRF.iv["Y","Y",],
                     Y.to.G.lower68.iv = IRF.iv["Y","G",]-IRF.se.iv[,"YG"],
                     Y.to.T.lower68.iv = IRF.iv["Y","T",]-IRF.se.iv[,"YT"],
                     Y.to.Y.lower68.iv = IRF.iv["Y","Y",]-IRF.se.iv[,"YY"],
                     Y.to.G.upper68.iv = IRF.iv["Y","G",]+IRF.se.iv[,"YG"],
                     Y.to.T.upper68.iv = IRF.iv["Y","T",]+IRF.se.iv[,"YT"],
                     Y.to.Y.upper68.iv = IRF.iv["Y","Y",]+IRF.se.iv[,"YY"]
                     )


figA2A.1 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=G.to.G),color="red")+
  geom_line(aes(y=G.to.G.iv),color="black")+
  geom_line(aes(y=G.to.G.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=G.to.G.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=G.to.G.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=G.to.G.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="g response", title="")+
  theme(legend.position = "none")+
  fig_thm

figA2A.2 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=T.to.G),color="red")+
  geom_line(aes(y=T.to.G.iv),color="black")+
  geom_line(aes(y=T.to.G.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=T.to.G.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=T.to.G.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=T.to.G.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="tax response", title="Goverment spending")+
  theme(legend.position = "none")+
  fig_thm

figA2A.3 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=Y.to.G),color="red")+
  geom_line(aes(y=Y.to.G.iv),color="black")+
  geom_line(aes(y=Y.to.G.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=Y.to.G.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=Y.to.G.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=Y.to.G.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="")+
  theme(legend.position = "none")+
  fig_thm

#Figure A2A
fig.A2A = plot_grid(figA2A.1,figA2A.2,figA2A.3,ncol=3)


figA2B.1 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=G.to.T),color="red")+
  geom_line(aes(y=G.to.T.iv),color="black")+
  geom_line(aes(y=G.to.T.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=G.to.T.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=G.to.T.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=G.to.T.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="g response", title="")+
  theme(legend.position = "none")+
  fig_thm

figA2B.2 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=T.to.T),color="red")+
  geom_line(aes(y=T.to.T.iv),color="black")+
  geom_line(aes(y=T.to.T.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=T.to.T.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=T.to.T.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=T.to.T.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="tax response", title="Tax shock")+
  theme(legend.position = "none")+
  fig_thm

figA2B.3 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=Y.to.T),color="red")+
  geom_line(aes(y=Y.to.T.iv),color="black")+
  geom_line(aes(y=Y.to.T.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=Y.to.T.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=Y.to.T.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=Y.to.T.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="")+
  theme(legend.position = "none")+
  fig_thm

#Figure A2B
fig.A2B = plot_grid(figA2B.1,figA2B.2,figA2B.3,ncol=3)

figA2C.1 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=G.to.Y),color="red")+
  geom_line(aes(y=G.to.Y.iv),color="black")+
  geom_line(aes(y=G.to.Y.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=G.to.Y.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=G.to.Y.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=G.to.Y.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="g response", title="")+
  theme(legend.position = "none")+
  fig_thm

figA2C.2 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=T.to.Y),color="red")+
  geom_line(aes(y=T.to.Y.iv),color="black")+
  geom_line(aes(y=T.to.Y.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=T.to.Y.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=T.to.Y.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=T.to.Y.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="tax response", title="Output shock")+
  theme(legend.position = "none")+
  fig_thm

figA2C.3 = ggplot(data=irfs.df,aes(x=horizon),group=1)+
  geom_line(aes(y=0))+
  geom_line(aes(y=Y.to.Y),color="red")+
  geom_line(aes(y=Y.to.Y.iv),color="black")+
  geom_line(aes(y=Y.to.Y.lower68),linetype="longdash",color="red")+
  geom_line(aes(y=Y.to.Y.upper68),linetype="longdash",color="red")+
  geom_line(aes(y=Y.to.Y.lower68.iv),linetype="longdash",color="black")+
  geom_line(aes(y=Y.to.Y.upper68.iv),linetype="longdash",color="black")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Quarters", y="y response", title="")+
  theme(legend.position = "none")+
  fig_thm

#Figure A2C
fig.A2C = plot_grid(figA2C.1,figA2C.2,figA2C.3,ncol=3)
