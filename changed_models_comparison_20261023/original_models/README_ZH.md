# 四个模型本来的样子

- base_layer_models.png / pdf：本轮实际使用的四层基础介质、真实Q、速度、密度和层厚。
- clean_reference_gathers.png / pdf：四个模型正常AVO的精确参考道集，未衰减、未加噪、未补偿，统一色标。
- Changed_*_all_scenes.png / pdf：每模型五种场景的精确纯净参考，模型内统一色标。
- layer_parameters.csv：完整分层参数；model_parameters.json：原实验参数计划。
- data/：20个模型/场景组合的参考数据，CSV和NPZ，均直接取自已经运行保存的truth字段；CSV第一列为秒，随后为11道。

基础介质为一维四层模型，界面深度由实验采用的 h=Δt×Vp/2 得到。第四层为半空间，地层图只截画至1300米；横向宽度用于示意。
断层、弱反射和尖灭场景在参考反射事件的时间或系数上作变化，再通过一维衰减算子产生观测；不是含断层或尖灭体的完整二维波动方程正演。
参考道集采用Shuey两项AVO与Ricker子波，已按反射时间对齐；不包含几何扩散、透射、多次波或NMO拉伸。
WeakSection第二反射系数乘0.08；FaultSection有8毫秒正弦起伏，第7至11道的深部两反射延迟40毫秒；PinchoutSection第二、第三反射异号并随横向靠拢，最后三道完全抵消。
薄层极性反转模型WeakSection中央道的第二、第三反射重合，不能分别恢复振幅。
