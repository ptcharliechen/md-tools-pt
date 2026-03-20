from numpy import arange
from collections.abc import Iterable
import matplotlib.pyplot as plt

class Plot:
    def __init__(self):
        self._data = []
        self._label = None
        self._fig = plt.figure()
        self._host = self._fig.add_subplot(111)
        self._label_flag, self._xlabel_flag, self._ylabel_flag = False, False, False
        self._xticks_flag, self._yticks_flag = False, False
        self._xlim_flag, self._ylim_flag = False, False
        self._xmin, self._xmax = None, None
        self._ymin, self._ymax = None, None
        plt.subplots_adjust(right=0.95, top=0.95, bottom=0.118)
    @property
    def data(self):
        return self._data
    @property
    def line_label(self):
        return self._label_flag if not self._label_flag else self._label
    @line_label.setter
    def line_label(self, label):
        assert len(label) == len(self._data), "Label length does not match data length."
        self._label = label
        self._label_flag = True
    @property
    def xlabel(self):
        return self._xlabel_flag
    @xlabel.setter
    def xlabel(self, xlabel):
        self._host.set_xlabel(xlabel, {"fontsize": 15})
        self._xlabel_flag = True
    @property
    def ylabel(self):
        return self._ylabel_flag
    @ylabel.setter
    def ylabel(self, ylabel):
        self._host.set_ylabel(ylabel, {"fontsize": 15})
        self._ylabel_flag = True
    @property
    def xticks(self):
        return self._xticks_flag
    @xticks.setter
    def xticks(self, xticks):
        data = []
        for iterator in range(len(self._data)):
            data += self._data[iterator][0]
        rng = max(data) - min(data)
        plt.xticks(arange(min(data), max(data)+rng/100, rng/(len(xticks)-1)), xticks, fontsize=12)
        self._xticks_flag = True
    @property
    def yticks(self):
        return self._yticks_flag
    @yticks.setter
    def yticks(self, yticks):
        data = []
        for iterator in range(len(self._data)):
            data += self._data[iterator][1]
        rng = max(data) - min(data)
        self._yticks_flag = True
        self._host.yticks(arange(min(data), max(data)+rng/100, rng/len(yticks)), yticks, fontsize=12)
    def auto_xticks(self):
        from numpy import arange, ceil
        x_len = self._data[0][0]
        if len(x_len) < 1001:
            xticks = None
        elif len(x_len) > 1000 and len(x_len) < 5001:
            xticks = list(arange(0, ceil(len(x_len/1000))))
        elif len(x_len) > 5000 and len(x_len) < 10001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 2))
        elif len(x_len) > 10000 and len(x_len) < 25001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 5))
        elif len(x_len) > 25000 and len(x_len) < 80001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 10))
        elif len(x_len) > 80000 and len(x_len) < 200001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 25))
        elif len(x_len) > 200000 and len(x_len) < 1000001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 100))
        elif len(x_len) > 1000000:
            xticks = list(arange(0, ceil(len(x_len)/1000000)))
        if xticks is not None:
            self.xticks(xticks)
    @property
    def xlim(self):
        return self._xlim_flag
    @xlim.setter
    def xlim(self, xlim):
        if not isinstance(xlim, Iterable):
            self._xmax = xlim
        elif xlim[0] is None and xlim[1] is not None:
            self._xmax = xlim[1]
        elif xlim[1] is None and xlim[0] is not None:
            self._xmin = xlim[0]
        elif xlim[0] is not None and xlim[1] is not None:
            self._xmin, self._xmax = xlim
        else:
            raise ValueError("Invalid value for xlim.")
        self._xlim_flag = True
    @property
    def ylim(self):
        return self._ylim_flag
    @ylim.setter
    def ylim(self, ylim):
        if not isinstance(ylim, Iterable):
            self._ymax = ylim
        elif ylim[0] is None and ylim[1] is not None:
            self._ymax = ylim[1]
        elif ylim[1] is None and ylim[0] is not None:
            self._ymin = ylim[0]
        elif ylim[0] is not None and ylim[1] is not None:
            self._ymin, self._ymax = ylim
        else:
            raise ValueError("Invalid value for ylim.")
        self._ylim_flag = True
    def append(self, x, y):
        self._data.append((x, y))
    def plot(self, legend=True, palette=None):
        from seaborn import color_palette
        if palette is None:
            palette = color_palette("Set2", 8)
            if len(self.data) > 8:
                palette += color_palette("husl", len(self.data)-8)
        if len(self._data) < 6:
            fontsize = 14
        elif len(self._data) < 9:
            fontsize = 12
        elif len(self._data) < 11:
            fontsize = 10
        else:
            fontsize = 8
        if legend:
            for idx, data in enumerate(self._data):
                if self._label is None:
                    self._host.plot(data[0], data[1], color=palette[idx], label=str(idx), linewidth=2.2)
                else:
                    self._host.plot(data[0], data[1], color=palette[idx], label=self._label[idx], linewidth=2.2)
            self._host.legend(loc="best", fontsize=fontsize)
        if not self._xticks_flag:
            self._host.tick_params(axis='x', labelsize=14)
        if not self._yticks_flag:
            self._host.tick_params(axis='y', labelsize=14)
        for idx, data in enumerate(self._data):
            self._host.plot(data[0], data[1], color=palette[idx], linewidth=2.2)
        tmp_x, tmp_y = [], []
        for i in range(len(self._data)):
            tmp_x += list(self._data[i][0])
            tmp_y += list(self._data[i][1])
        buffer = (max(tmp_y) - min(tmp_y)) * 0.05
        if self._xmin is None:
            self._xmin = min(tmp_x)
        if self._xmax is None:
            self._xmax = max(tmp_x)
        if self._ymin is None:
            self._ymin = min(tmp_y)-buffer
        if self._ymax is None:
            self._ymax = max(tmp_y)+buffer
        self._host.set_xlim(self._xmin, self._xmax)
        self._host.set_ylim(self._ymin, self._ymax)
    def save(self, name):
        from kit.fundamental import Args
        from os import getcwd
        if ".png" not in name:
            name += ".png"
        Args.same_name(getcwd(), name)
        self._fig.savefig(name, dpi=600, format="png")

class TwinPlot(Plot):
    def __init__(self):
        super().__init__()
        self._left_data, self._right_data = [], []
        self._fig = plt.figure()
        self._host = self._fig.add_subplot(111)
        self._twin = self._host.twinx()
        self._left_ylim_flag, self._right_ylim_flag = False, False
        self._left_ymin, self._left_ymax = None, None
        self._right_ymin, self._right_ymax = None, None
        plt.subplots_adjust(top=0.95, bottom=0.118)
    @property
    def data(self):
        return self._data
    @property
    def left_ylabel(self):
        return self._ylabel_flag % 2
    @left_ylabel.setter
    def left_ylabel(self, ylabel):
        self._ylabel_flag += 2
        self._host.set_ylabel(ylabel, {"fontsize": 15})
    @property
    def right_ylabel(self):
        return self._ylabel_flag // 2
    @right_ylabel.setter
    def right_ylabel(self, ylabel):
        self._ylabel_flag += 1
        self._twin.set_ylabel(ylabel, {"fontsize": 15})
    @property
    def left_yticks(self):
        return self._yticks_flag % 2
    @left_yticks.setter
    def left_yticks(self, yticks):
        data = []
        for iterator in range(len(self._left_data)):
            data += self._data[iterator][1]
        rng = max(data) - min(data)
        self._yticks_flag += 2
        self._host.yticks(arange(min(data), max(data)+rng/100, rng/len(yticks)), yticks, fontsize=14)
    @property
    def right_yticks(self):
        return self._yticks_flag // 2
    @right_yticks.setter
    def right_yticks(self, yticks):
        data = []
        for iterator in range(len(self._right_data)):
            data += self._data[iterator][1]
        rng = max(data) - min(data)
        self._yticks_flag += 1
        self._twin.yticks(arange(min(data), max(data)+rng/100, rng/len(yticks)), yticks, fontsize=14)
    @property
    def left_ylim(self):
        return self._left_ylim_flag
    @left_ylim.setter
    def left_ylim(self, ylim):
        if not isinstance(ylim, Iterable):
            self._left_ymax = ylim
        elif ylim[0] is None and ylim[1] is not None:
            self._left_ymax = ylim[1]
        elif ylim[1] is None and ylim[0] is not None:
            self._left_ymin = ylim[0]
        elif ylim[0] is not None and ylim[1] is not None:
            self._left_ymin, self._left_ymax = ylim
        else:
            raise ValueError("Invalid value for ylim.")
        self._left_ylim_flag = True
    @property
    def right_ylim(self):
        return self._right_ylim_flag
    @right_ylim.setter
    def right_ylim(self, ylim):
        if not isinstance(ylim, Iterable):
            self._right_ymax = ylim
        elif ylim[0] is None and ylim[1] is not None:
            self._right_ymax = ylim[1]
        elif ylim[1] is None and ylim[0] is not None:
            self._right_ymin = ylim[0]
        elif ylim[0] is not None and ylim[1] is not None:
            self._right_ymin, self._right_ymax = ylim
        else:
            raise ValueError("Invalid value for ylim.")
        self._right_ylim_flag = True
    def append(self, x, y, axis="left"):
        if axis == "left":
            self._left_data.append((x, y))
            self._data += self._left_data
        elif axis == "right":
            self._right_data.append((x, y))
            self._data += self._right_data
    def plot(self, palette=None, legend=True):
        from seaborn import color_palette
        if palette is None:
            palette = color_palette("Set2", 8)
            if len(self.data) > 8:
                palette += color_palette("husl", len(self.data)-8)
        if len(self._left_data) < 6:
            left_fontsize = 14
        elif len(self._left_data) < 9:
            left_fontsize = 12
        elif len(self._left_data) < 11:
            left_fontsize = 10
        else:
            left_fontsize = 8
        if len(self._right_data) < 6:
            right_fontsize = 14
        elif len(self._right_data) < 9:
            right_fontsize = 12
        elif len(self._right_data) < 11:
            right_fontsize = 10
        else:
            right_fontsize = 8
        if not self._xticks_flag:
            self._host.tick_params(axis='x', labelsize=14)
            self._twin.tick_params(axis='x', labelsize=14)
        if not self._yticks_flag:
            self._host.tick_params(axis='y', labelsize=14)
            self._twin.tick_params(axis='y', labelsize=14)
        if legend:
            for idx, data in enumerate(self._left_data):
                if len(self._label) == 0:
                    self._host.plot(data[0], data[1], color=palette[idx], label=str(idx), linewidth=2.2)
                else:
                    self._host.plot(data[0], data[1], color=palette[idx], label=self._label[idx], linewidth=2.2)
            for idx, data in enumerate(self._right_data, start=len(self._left_data)):
                if len(self._label) == 0:
                    self._twin.plot(data[0], data[1], color=palette[idx], label=str(idx), linewidth=2.2)
                else:
                    self._twin.plot(data[0], data[1], color=palette[idx], label=self._label[idx-len(self._left_data)+1], linewidth=2.2)
            self._host.legend(loc="upper left", fontsize=left_fontsize)
            self._twin.legend(loc="upper right", fontsize=right_fontsize)
        else:
            for idx, data in enumerate(self._left_data):
                self._host.plot(data[0], data[1], color=palette[idx], linewidth=2.2)
            for idx, data in enumerate(self._right_data, start=len(self._left_data)):
                self._twin.plot(data[0], data[1], color=palette[idx], linewidth=2.2)
        tmp_x, tmp_y_left, tmp_y_right = [], [], []
        for i in range(len(self._left_data)):
            tmp_x += list(self._data[i][0])
            tmp_y_left += list(self._left_data[i][1])
            tmp_y_right += list(self._right_data[i][1])
        left_buffer, right_buffer = (max(tmp_y_left) - min(tmp_y_left)) * 0.05, (max(tmp_y_right) - min(tmp_y_right)) * 0.05
        if self._xmin is None:
            self._xmin = min(tmp_x)
        if self._xmax is None:
            self._xmax = max(tmp_x)
        self._host.set_xlim(self._xmin, self._xmax)
        self._twin.set_xlim(self._xmin, self._xmax)
        if self._left_ymin is None:
            self._left_ymin = min(tmp_y_left)-left_buffer
        if self._left_ymax is None:
            self._left_ymax = max(tmp_y_left)+left_buffer
        if self._right_ymin is None:
            self._right_ymin = min(tmp_y_right)-right_buffer
        if self._right_ymax is None:
            self._right_ymax = max(tmp_y_right)+right_buffer
        self._host.set_ylim(self._left_ymin, self._left_ymax)
        self._twin.set_ylim(self._right_ymin, self._right_ymax)

def plot(fig, file_name, ticks=False, line_label=True, xlabel=True, ylabel=True, palette=None):
    if ticks:
        from numpy import arange, ceil
        while(1):
            tmp = input("Adjust X axis ticks, or auto mode (y/n/a) [a]: ")
            if tmp == '' or tmp.lower() == 'a':
                x_len = fig.data[0][0]
                if len(x_len) > 1000 and len(x_len) < 5001:
                    fig.xticks = list(arange(0, ceil(len(x_len/1000))))
                elif len(x_len) > 5000 and len(x_len) < 10001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 2))
                elif len(x_len) > 10000 and len(x_len) < 25001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 5))
                elif len(x_len) > 25000 and len(x_len) < 80001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 10))
                elif len(x_len) > 80000 and len(x_len) < 200001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 25))
                elif len(x_len) > 200000 and len(x_len) < 1000001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 100))
                elif len(x_len) > 1000000:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000000)))
                break
            elif tmp.lower() == 'n':
                break
            else:
                print("Warning: Input error.")
    if not line_label:
        name_flag = False
    if line_label and not fig.line_label:
        while(1):
            tmp = input("Name the lines (y/n) [n]: ")
            if tmp == '':
                name_flag = False
                break
            elif tmp.lower() == 'y':
                name_flag = True
                break
            elif tmp.lower() == 'n':
                name_flag = False
                break
            else:
                print("Warning: Input error.")
        if name_flag:
            label = []
            if isinstance(fig, Plot):
                for i in range(len(fig.data)):
                    line = input(f"Line {i+1}: ")
                    label.append(line)
            elif isinstance(fig, TwinPlot):
                for i in range(len(fig.data[0])):
                    line = input(f"Left line {i+1}: ")
                    label.append(line)
                    line = input(f"Right line {i+1}: ")
                    label.append(line)
            fig.line_label = label
    elif fig.line_label:
        name_flag = True
    if xlabel and not fig.xlabel:
        while(1):
            tmp = input("Name X axis (y/n) [n]: ")
            if tmp == '':
                break
            elif tmp.lower() == 'y':
                fig.xlabel = input("X axis: ")
                break
            elif tmp.lower() == 'n':
                break
            else:
                print("Warning: Input error.")
    if isinstance(fig, Plot):
        if not fig.ylabel:
            while(1):
                tmp = input("Name Y axis (y/n) [n]: ")
                if tmp == '':
                    break
                elif tmp.lower() == 'y':
                    fig.ylabel = input("Y axis: ")
                    break
                elif tmp.lower() == 'n':
                    break
                else:
                    print("Warning: Input error.")
    elif ylabel and isinstance(fig, TwinPlot):
        if not fig.left_ylabel:
            while(1):
                tmp = input("Name left Y axis (y/n) [n]: ")
                if tmp == '':
                    break
                elif tmp.lower() == 'y':
                    fig.left_ylabel = input("Left Y axis: ")
                    break
                elif tmp.lower() == 'n':
                    break
                else:
                    print("Warning: Input error.")
        elif ylabel and not fig.right_ylabel:
            while(1):
                tmp = input("Name right Y axis (y/n) [n]: ")
                if tmp == '':
                    break
                elif tmp.lower() == 'y':
                    fig.right_ylabel = input("Right Y axis: ")
                    break
                elif tmp.lower() == 'n':
                    break
                else:
                    print("Warning: Input error.")
    if palette is None:
        while(1):
            tmp = input("Adjust color (y/n) [n]: ")
            if tmp == '':
                palette = None
                break
            elif tmp.lower() == 'y':
                palette = []
                for i in range(len(fig.data)):
                    tmp = input(f"Color {i+1}: ")
                    palette.append(tmp)
                break
            elif tmp.lower() == 'n':
                palette = None
                break
            else:
                print("Warning: Input error.")
    fig.plot(palette=palette, legend=name_flag)
    fig.save(file_name)
