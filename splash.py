image = """                                             
\t\t\t\t\t                      `-:+--                        
\t\t\t\t\t                  -sdNMMNhsMm+                      
\t\t\t\t\t                `syyhs+-` `oNMm-                    
\t\t\t\t\t               .dMN-        `oNN/                   
\t\t\t\t\t               ydo.        `::+dy                   
\t\t\t\t\t              :MMy         +/oyh/                   
\t\t\t\t\t              sMMs         `                        
\t\t\t\t\t              ohh:                                  
\t\t\t\t\t              +NMN.     :/.oys::                    
\t\t\t\t\t              .mhyo`   `my.::..::                   
\t\t\t\t\t               .hdds-``yd` /d/hh+/-                 
\t\t\t\t\t                -mNMNmsy. .m/ .-``:.                
\t\t\t\t\t         `s+/:::-MMMMMMMd/so  od:yy+/:              
\t\t\t\t\t        :yo:/+++-NMMyhhddho. /m-``.``::             
\t\t\t\t\t       oNy:/.`  `hNyhMMMMMMm/h- oy:sso::`           
\t\t\t\t\t      :o::ooossss-ymMMMshddmh- /m:`:+/-:+           
\t\t\t\t\t      ../m+ ``    oMMNymMNNNNN:d- .ymmdo.           
\t\t\t\t\t       `y- .hsssss+odymMMhdmmms. :d+..:hm-`         
\t\t\t\t\t        + +d+`````` +MMMNomdmmm/-mh`    hyyo`       
\t\t\t\t\t        `:h/  .hsosso+NmsNMMMMMooo`     .oMN/       
\t\t\t\t\t         o- `/y+....` -+NMMMMMMNoss`     -hyys+-    
\t\t\t\t\t         `. od/         :osmMNy: .+/      oMMMMNo   
\t\t\t\t\t           `y`     `  `:+hooy.`    .     `MMhMMMM+  
\t\t\t\t\t            :     -+h+oy+- om:`          `NM sMMMs  
\t\t\t\t\t                  `dMm      .:.           :m-`MMM:  
\t\t\t\t\t                   .dNo``                  `.`MN+   
\t\t\t\t\t                     /ydhs/o:`               sh-    
\t\t\t\t\t                       `+dMMd++syy+`        `.      
\t\t\t\t\t                         -//hMMNyymm-               
\t\t\t\t\t                           `NMMMy. -o               
\t\t\t\t\t                            -dMMMMds+//:            
\t\t\t\t\t                              `:/ooso+:`            
                                                  
"""

pretitle = """
\t\t____________________________________________________________________________________________________________
"""

title = """
\t\t___/\/\/\/\/\____/\/\/\/\/\____/\/\/\/\____/\/\/\/\/\____/\/\/\/\/\____/\/\/\/\____/\/\/\/\____/\/\____/\/\_
\t\t_/\/\__________/\/\__________/\/\____/\/\__/\/\____/\/\__/\/\____/\/\____/\/\____/\/\____/\/\__/\/\/\__/\/\_
\t\t___/\/\/\/\____/\/\__________/\/\____/\/\__/\/\/\/\/\____/\/\/\/\/\______/\/\____/\/\____/\/\__/\/\/\/\/\/\_
\t\t_________/\/\__/\/\__________/\/\____/\/\__/\/\__/\/\____/\/\____________/\/\____/\/\____/\/\__/\/\__/\/\/\_
\t\t_/\/\/\/\/\______/\/\/\/\/\____/\/\/\/\____/\/\____/\/\__/\/\__________/\/\/\/\____/\/\/\/\____/\/\____/\/\_
\t\t____________________________________________________________________________________________________________
"""

print "\n"
#print image
#print title
CSI = "\x1B[49m\x1B["
print CSI+"92m"+pretitle+CSI+"0m"
print CSI+"93m"+image+CSI+"0m"
print CSI+"92m"+title+CSI+"0m"
print "\n"
