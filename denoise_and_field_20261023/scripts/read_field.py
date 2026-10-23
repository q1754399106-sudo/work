from pathlib import Path
import struct,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'field_data'
data=(OUT/'wz.25.H').read_bytes();offset=data.find(b'\x0c\x0c\x04')+3;assert offset==1284
header=data[:offset-3].decode('ascii');(OUT/'wz25_header.txt').write_text(header,encoding='utf-8')
a=np.frombuffer(data[offset:],dtype='>f4').astype(float).reshape(97,2000).T;assert np.isfinite(a).all()
x=(np.arange(97)-48)*50.;t=.002+np.arange(2000)*.002;np.savez_compressed(OUT/'wz25_raw.npz',data=a,t=t,x=x,dt=.002)
meta=dict(wz25=dict(area='Alberta, Canada',source='Dynamite',sample_count=2000,trace_count=97,dt_s=.002,t0_s=.002,offset_spacing_m=50,offset_range_m=[-2400,2400],center_trace_note='Header says a zero-valued inner-offset trace was inserted at trace49; retained, detected as dead for processing',geometry_note='Historical header notes mention uncertainty in spacing, retained declared50m; not independently surveyed',header_file='wz25_header.txt',source_sha256=hashlib.sha256(data).hexdigest()))
data=(OUT/'L23535.SGY').read_bytes();header=data[:3200].decode('cp500');(OUT/'alaska_text_header.txt').write_text('\n'.join(header[j:j+80] for j in range(0,3200,80)),encoding='utf-8')
dt=struct.unpack_from('>H',data,3216)[0]*1e-6;ns=struct.unpack_from('>H',data,3220)[0];fmt=struct.unpack_from('>H',data,3224)[0];assert dt==.002 and ns==3000 and fmt==1
step=240+ns*4;assert (len(data)-3600)%step==0;n=(len(data)-3600)//step
groups={}
for j in range(n):
 pos=3600+j*step;fldr=struct.unpack_from('>i',data,pos+8)[0];chan=struct.unpack_from('>i',data,pos+12)[0]
 if fldr not in [111,120,130] or chan>96:continue
 words=np.frombuffer(data,dtype='>u4',count=ns,offset=pos+240).astype(np.uint32)
 sign=np.where(words>>31,-1.,1.);exponent=((words>>24)&0x7f).astype(int)-64;fraction=(words&0xffffff)/16777216.
 samples=sign*fraction*np.power(16.,exponent);samples[words==0]=0;assert np.isfinite(samples).all();groups.setdefault(fldr,[]).append((chan,samples))
for shot,gather in groups.items():
 gather.sort();assert [k for k,v in gather]==list(range(1,97));a=np.column_stack([v for k,v in gather]);x=(-5225+110*np.arange(96))*.3048;t=np.arange(ns)*dt
 np.savez_compressed(OUT/f'alaska_{shot}_raw.npz',data=a,t=t,x=x,dt=dt)
meta['alaska']=dict(origin='USGS NPRA Alaska, line31-81,1981, raw demultiplexed SEG-Y',format='Big endian IBM32',dt_s=dt,ns=ns,file_trace_count=n,field_record_range=[101,130],selected_records=[111,120,130],physical_channels=96,aux_channels_excluded=[97,98,99,100,101],offset_spacing_m=33.528,offset_range_m=[-1592.58,1592.58],geometry_note='Reconstructed from upstream Madagascar processing script o2=-5225ft,d2=110ft; raw offset header is zero. Not navigation coordinates.',source_sha256=hashlib.sha256(data).hexdigest())
(OUT/'parsed_metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(meta,ensure_ascii=True,indent=2))
