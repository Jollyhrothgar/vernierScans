select
  trigger.runnumber,run.triggerconfig,trigger.name,trigger.bitnb,trigger.scaleraerraw,trigger.scaleraerlive
from
  run,trigger
where
  run.runnumber = trigger.runnumber
  and (trigger.bitnb = 11 or trigger.bitnb = 1 or trigger.bitnb = 2 or trigger.bitnb = 0)
  and (run.runnumber = 359711 
	or run.runnumber = 365866 
	or run.runnumber = 362492 
	or run.runnumber = 360879 
	or run.runnumber = 366604
	or run.runnumber = 362138
  )
