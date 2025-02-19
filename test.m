scene = sceneCreate;    scene = sceneSet(scene,'fov',10);
oi = oiCreate('wvf');   oi = oiCompute(oi,scene,'crop',true);
sensor = sensorCreate('imx363');  sensor = sensorSet(sensor,'fov',10,oi);
sensor = sensorCompute(sensor,oi);
sensorWindow(sensor);
