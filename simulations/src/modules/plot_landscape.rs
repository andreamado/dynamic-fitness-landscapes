use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufWriter, Write},
    process::Command
};

use super::genotype::Genotype;

#[derive(Clone)]
pub enum Color {
    RGB(i32, i32, i32),
    Hex(String)
}

impl Color {
    fn to_rgb(&self) -> Self {
        match self {
            Color::RGB(..) => (*self).clone(),
            Color::Hex(s)  => {
                let r = i32::from_str_radix(&s[1..3], 16).unwrap();
                let g = i32::from_str_radix(&s[3..5], 16).unwrap();
                let b = i32::from_str_radix(&s[5..7], 16).unwrap();
                Color::RGB(r,g,b)
            }
        }
    }

    fn to_hex(&self) -> Self {
        match self {
            Color::RGB(r,g,b) => {
                Color::Hex(format!("#{:X}{:X}{:X}", r, g, b))
            },
            Color::Hex(_) => (*self).clone()
        }
    }

    fn as_string(&self) -> String {
        match self.to_hex() {
            Color::Hex(s) => s.clone(),
            _ => unreachable!()
        }
    }

    fn as_tuple(&self) -> (i32, i32, i32) {
        match self.to_rgb() {
            Color::RGB(r,g,b) => (r,g,b),
            _ => unreachable!()
        }
    }

    fn gradient(c1: &Color, c2: &Color, pos: f64) -> Color {
        let dif = (*c2).clone() - (*c1).clone();
        (*c1).clone() + dif * pos
    }
}

impl std::ops::Add<Color> for Color {
    type Output = Color;

    fn add(self, rhs: Color) -> Color {
        let a = self.as_tuple();
        let b = rhs.as_tuple();
        Self::RGB(a.0+b.0, a.1+b.1, a.2+b.2)
    }
}
impl std::ops::Sub<Color> for Color {
    type Output = Color;

    fn sub(self, rhs: Color) -> Color {
        let a = self.as_tuple();
        let b = rhs.as_tuple();
        Self::RGB(a.0-b.0, a.1-b.1, a.2-b.2)
    }
}
impl std::ops::Mul<f64> for Color {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        let a = self.as_tuple();
        Color::RGB((a.0 as f64 * rhs) as i32, (a.1 as f64 * rhs) as i32, (a.2 as f64 * rhs) as i32)
    }
}

fn get_connections<const L: usize>(genotypes: &Vec<Genotype<L>>) -> Vec<(&Genotype<L>, &Genotype<L>)> {
    let mut connections = Vec::with_capacity(genotypes.len()*genotypes.len());
    for g1 in genotypes {
        for g2 in genotypes {
            if g1.sum() == g2.sum() + 1 && g1.n_differences(g2) == 1 {
                connections.push((g1, g2));
            }
        }
    }
    connections.sort();
    connections.dedup();
    connections
}

fn line(beg: (f64, f64), end: (f64, f64), color: &str, thickness: f64, opacity: f64) -> String {
    format!(
r##"    <line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" style="stroke: {color}; stroke-width:{thickness}%%; stroke-opacity={opacity};" />
"##,
        x1=beg.0, y1=beg.1, x2=end.0, y2=end.1, color=color, thickness=thickness, opacity=opacity
    )
}

fn text(text: &str, pos: (f64, f64), font_size: f64, align: &str) -> String {
    format!(
r##"    <text x="{x:.1}" y="{y:.1}" text-anchor="{align}" style="font-size:{font_size}pt;">{text}</text>
"##,
        text=text, x=pos.0, y=pos.1, font_size=font_size, align=align
    )
}

fn rectangle(pos: (f64, f64), size: (f64, f64), color: &str, opacity: f64, round_corners: f64) -> String {
    format!(
r#"    <rect x="{x:.1}" y="{y:.1}" rx="{r:.1}" ry="{r:.1}" width="{width:.1}" height="{height:.1}" fill="{background}" fill-opacity="{opacity}" />
"#,
        x=pos.0, y=pos.1, r=round_corners,
        width=size.0, height=size.1,
        background=color, opacity=opacity
    )
}

pub enum Ticks {
    Number(usize),
    List(Vec<f64>),
    LabeledList(Vec<(f64, String)>)
}


pub struct FitnessLandscapePlot<'a, const L: usize> {
    landscape: &'a HashMap::<Genotype<L>, f64>,
    landscape_std: Option<&'a HashMap::<Genotype<L>, f64>>,
    colors: Option<&'a HashMap::<Genotype<L>, f64>>,
    genotypes: Vec<Genotype<L>>,

    pub size: (f64, f64),
    pub margins: ((f64, f64), (f64, f64)),
    pub ylims: (f64, f64),

    pub ticks: Ticks,
    pub tick_size: f64,
    pub tick_distance: f64,
    pub tick_font_size: f64,
    pub tick_precision: usize,

    pub marker_shift: f64,
    pub marker_width: f64,
    pub marker_height: f64,
    pub marker_radius: f64,

    pub marker_color: [Color; 2],
    pub connection_colors: [&'a str; 2],
    pub marker_gene_color: [&'a str; 2],
    pub axis_color: &'a str,
    pub background_color: &'a str,

    pub axis_tickness: f64,

    pub connections: bool,

    pub render: &'a str,

    pub labels_bottom: bool
}

impl<'a, const L: usize> FitnessLandscapePlot<'a, L> {
    pub fn new(landscape: &'a HashMap::<Genotype<L>, f64>, landscape_std: Option<&'a HashMap::<Genotype<L>, f64>>, colors: Option<&'a HashMap::<Genotype<L>, f64>>) -> Self {
        let mut genotypes: Vec<Genotype<L>> = landscape.keys().map(|g| (*g).clone()).collect();
        genotypes.sort_by(|g1, g2|
            g1.order().partial_cmp(&g2.order()).unwrap()
        );

        let flmax = match landscape_std {
            Some(std) => landscape.iter().map(|(g, f)| f + std[g]).max_by(
                            |fa, fb| fa.partial_cmp(fb).expect("Tried to compare a NaN")
                         ).unwrap(),
            None      => *landscape.values()
                                   .max_by(
                                       |fa, fb| fa.partial_cmp(fb).expect("Tried to compare a NaN")
                                   ).unwrap()
        };
        let flmin = match landscape_std {
            Some(std) => landscape.iter().map(|(g, f)| f - std[g]).max_by(
                            |fa, fb| fb.partial_cmp(fa).expect("Tried to compare a NaN")
                         ).unwrap(),
            None      => *landscape.values()
                                   .max_by(
                                       |fa, fb| fb.partial_cmp(fa).expect("Tried to compare a NaN")
                                   ).unwrap()
        };

        let r = match colors {
            Some(_) => 400.,
            None    => 200.
        };

        Self {
            landscape,
            landscape_std,
            genotypes,
            colors,

            size: (2000., 1200.),
            margins: ((200., r), (120., 50.)),
            ylims: (flmin, flmax),

            ticks: Ticks::Number(4),
            tick_size: 10.,
            tick_distance: 10.,
            tick_font_size: 28.,
            tick_precision: 2,

            marker_shift: 10.,
            marker_width: 10.,
            marker_height: 30.,
            marker_radius: 6.,

            marker_color: [Color::Hex("#AAAAAA".to_string()), Color::Hex("#DC143C".to_string())],

            axis_color: "black",
            background_color: "white",

            axis_tickness: 1.,

            connections: true,
            connection_colors: ["#FFB3BF", "#CCCCFF"],
            marker_gene_color: ["#B2B2B2", "#0A66C2"],

            render: "",

            labels_bottom: true
        }
    }

    pub fn plot(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        // brew install librsvg
        // https://superuser.com/questions/134679/command-line-application-for-converting-svg-to-png-on-mac-os-x

        let ((l, r), (b, t)) = self.margins;
        let (w, h) = self.size;

        let mut graph = String::new();

        ///////////////////////////////////////////////////////////////////////////////////////////
        // header
        graph.push_str(format!(
r#"<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" style="background-color:{background_color}">
"#,
        width=w, height=h, background_color=self.background_color).as_str());
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // axes
        graph.push_str("    <!-- Fill the background -->\n");
        graph.push_str(r#"    <rect width="100%" height="100%" fill="white" class="background"/>"#);
        graph.push_str("\n");
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // axes
        graph.push_str("\n    <!-- Draw the axes -->\n");
        // left axis
        let (beg, end) = ((l, t), (l, h - b));
        graph.push_str(line(beg, end, self.axis_color, self.axis_tickness, 1.).as_str());

        // right axis
        let (beg, end) = ((w - r, t), (w - r, h - b));
        graph.push_str(line(beg, end, self.axis_color, self.axis_tickness, 1.).as_str());

        // horizontal axis
        let (beg, end) = ((l, self.to_y(1.)), (w - r, self.to_y(1.)));
        graph.push_str(format!(
r##"    <line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="black" stroke-width="1" stroke-opacity="1" stroke-dasharray="4" class="xaxis"/>
"##,
        x1=beg.0, y1=beg.1, x2=end.0, y2=end.1
        ).as_str());

        //ffmpeg -i video/%04d.svg -width 600 -vf format=yuv420p output.mp4

        graph.push_str(format!(
r##"    <text x="{x:.1}" y="{y:.1}" text-anchor="middle" transform="rotate(270 {x} {y})" style="font-size:{font_size}pt;" class="ylabel">Fitness</text>
"##,
        y = (h - t - b)/2. + t, x = l / 2. - self.tick_font_size*1.5, font_size = self.tick_font_size*1.5
        ).as_str());
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // ticks
        graph.push_str("\n");
        graph.push_str("    <!-- Draw the ticks -->\n");
        for (tick_pos, tick_label, tick_value) in self.generate_ticks().iter() {
            let tick_pos = *tick_pos;

            // left ticks
            let beg = (l,                  tick_pos);
            let end = (l + self.tick_size, tick_pos);
            graph.push_str(line(beg, end, self.axis_color, 1., 1.).as_str());

            let text_pos = (beg.0 - self.tick_distance, beg.1 + self.tick_font_size*0.3);
            graph.push_str(text(tick_label.as_str(), text_pos, self.tick_font_size, "end").as_str());

            // right ticks
            let beg = (w - r,                  tick_pos);
            let end = (w - r - self.tick_size, tick_pos);
            graph.push_str(line(beg, end, self.axis_color, 1., 1.).as_str());

            let mut x0 = beg.0 + self.tick_distance;
            if *tick_value > 0. { x0 += 3. };

            let text_pos = (x0, beg.1 + self.tick_font_size*0.3);
            graph.push_str(text(tick_label.as_str(), text_pos, self.tick_font_size, "start").as_str());
        }
        graph.push_str("\n");
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // connections
        if self.connections {
            let connections = get_connections(&self.genotypes);
            graph.push_str("    <!-- Draw the connections -->\n");

            if self.labels_bottom {
                let x_positions = self.generate_x_positions();

                for (i, s1) in self.genotypes.iter().enumerate() {
                    let f1 = self.landscape[s1];
                    let beg = (x_positions[i], self.to_y(f1));
                    for j in 0..L {
                        let mut s2 = (*s1).clone();
                        s2.mutate(j);
                        if s1.n_differences(&s2) == 1 && s2.sum() > s1.sum() {
                            let k = self.genotypes.iter().position(|g| g == &s2).unwrap();
                            let f2 = self.landscape[&s2];
                            let end = (x_positions[k], self.to_y(f2));
                            graph.push_str(line(beg, end, self.connection_colors[if f1 > f2 {0} else {1}], 0.1, 1.).as_str());
                        }
                    }
                }
            } else {
                for (seq1, seq2) in connections {
                    let f1 = self.landscape[seq1];
                    let f2 = self.landscape[seq2];

                    let beg = (
                        -self.marker_shift + l + (seq1.sum() as f64 + 0.25) * (w - l - r)/(L as f64 + 1.) - self.marker_width/2.,
                        self.to_y(f1)
                    );
                    let end = (
                         self.marker_shift + l + (seq2.sum() as f64 + 0.25) * (w - l - r)/(L as f64 + 1.) - self.marker_width/2. + self.marker_width*(L as f64),
                         self.to_y(f2)
                    );

                    graph.push_str(line(beg, end, self.connection_colors[if f1 > f2 {1} else {0}], 0.1, 1.).as_str());
                }
            }
            graph.push_str("\n");
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Genotype labels
        if self.labels_bottom {
            let block_positions = self.generate_block_positions();
            for x in block_positions.iter() {
                let pos  = (x.0, t);
                let size = (x.1-x.0, h - b - t + (L + 2) as f64 * self.marker_radius*2.5);
                graph.push_str(rectangle(pos, size, self.marker_color[0].as_string().as_str(), 0.2, self.marker_radius).as_str());
            }
            graph.push_str("\n");

            graph.push_str("    <!-- Draw the genotypes & sight guides -->\n");
            let x_positions = self.generate_x_positions();

            for (i, seq) in self.genotypes.iter().enumerate() {
                let x = x_positions[i];

                let yb = h - b + (L + 2) as f64 * (self.marker_radius*2.5);
                graph.push_str(format!(
r##"    <line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="black" stroke-width="0.5" stroke-opacity="0.2" stroke-dasharray="4" class="sight_guide"/>
"##, x1=x, y1=yb, x2=x, y2=t
                ).as_str());

                for (j, s) in seq.iter().enumerate() {
                    graph.push_str(format!(
r#"    <circle cx="{cx:.2}" cy="{cy:.2}" r="{r:.2}" fill="{color}" class="genotype_label" />
"#,
                        cx=x, cy=(h - b + (j + 2) as f64 * (self.marker_radius*2.5)), r=self.marker_radius,
                        color=self.marker_gene_color[*s as usize]
                    ).as_str());
                }
                graph.push_str("\n");
            }
            graph.push_str("    <!-- Draw fitness markers -->\n");
            for (i, x) in x_positions.iter().enumerate() {
                let x = *x;

                let g = &self.genotypes[i];
                let f = self.landscape[g];

                let mut occupation = 0.;
                let color = match self.colors {
                    Some(color_map) => {
                        occupation = *color_map.get(g).unwrap_or(&0.);
                        Color::gradient(&self.marker_color[0], &self.marker_color[1], occupation).as_string()
                    },
                    None => self.marker_color[0_usize].as_string()
                };

                let (pos, size) = match self.landscape_std {
                    Some(std) => {
                        let std = std[g];
                        ((x-self.marker_width/2., self.to_y(f + std)),
                         (self.marker_width, self.convert_height(2.*std)))
                    },
                    None => {
                        ((x, self.to_y(f)),
                         (0., 0.))
                    }
                };

                match self.landscape_std {
                    Some(_) => {
                        graph.push_str(rectangle(pos, size, color.as_str(), 1., self.marker_radius).as_str());
                    },
                    None => {
                        graph.push_str(format!(
r#"    <circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}" class="fitness_marker" />
"#,
                        cx=pos.0, cy=pos.1, r=self.marker_radius*(1.+occupation)*1.5, color=color).as_str());
                    }
                }
            }
        }
        graph.push_str("\n");
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // color label
        match self.colors {
            Some(_) => {
                graph.push_str("    <!-- Draw legend -->\n");
                let top  = h - 3.*b;
                let left = w - r*0.7;
                let right = w - r*0.05;
                let delta = right - left;

                graph.push_str(format!(
r##"    <text x="{x:.1}" y="{y:.1}" text-anchor="start" style="font-size:{font_size}pt;" class="legend_title">Frequency</text>
"##,
                y = top, x = left, font_size = self.tick_font_size*1.5
                ).as_str());

                let pos = (left - 1., top + self.tick_font_size*1.5);
                graph.push_str(rectangle(pos, (delta+1., 20.), self.marker_color[0].as_string().as_str(), 1., 0.).as_str());

                for i in 0..100 {
                    let x = i as f64 / 100.;
                    let pos = (left + x * delta, top + self.tick_font_size*1.5);
                    let size = (delta / 100.+1., 20.);
                    graph.push_str(rectangle(pos, size, Color::gradient(&self.marker_color[0], &self.marker_color[1], x).as_string().as_str(), 1., 0.).as_str());
                }
                graph.push_str("\n");

                graph.push_str("    <!-- Draw legend ticks -->\n");
                let n_ticks = 6;
                for i in 0..n_ticks {
                    let v = i as f64 / (n_ticks as f64 - 1.);
                    let x = left + v * delta;
                    let y = top + self.tick_font_size*2.5 + 20.;
                    graph.push_str(format!(
r##"    <text x="{x:.1}" y="{y:.1}" text-anchor="middle" style="font-size:{font_size}pt;" class="legend_tick">{v}</text>
"##,
                    y = y, x = x, font_size = self.tick_font_size*0.8, v = v
                    ).as_str());

                    let beg = (x, y      - self.tick_font_size);
                    let end = (x, y - 5. - self.tick_font_size);
                    graph.push_str(line(beg, end, "black", 1., 1.).as_str());
                }
                let beg = (left,  top + 20. + self.tick_font_size*1.5);
                let end = (right, top + 20. + self.tick_font_size*1.5);
                graph.push_str(line(beg, end, "black", 1., 1.).as_str());
            },
            None => {}
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////
        // footer
        graph.push_str("</svg>");
        ///////////////////////////////////////////////////////////////////////////////////////////

        let file = File::create(filename)?;
        let mut file = BufWriter::new(file);
        file.write(graph.as_bytes())?;
        file.flush()?;

        match self.render {
            "pdf" => {
                match Command::new("sh")
                    .args(["-c", format!("rsvg-convert -f pdf {} -o {}.pdf", filename, &filename[..(filename.len()-4)]).as_str(),])
                    .spawn() {
                      Ok(_) => {},
                      Err(_) => println!("Unable to generate png file (only svg generated). Is 'rsvg-convert' installed?")
                    }
            },
            "png" => {
                match Command::new("sh") 
                    .args(["-c", format!("rsvg-convert -f png {} -o {}.png", filename, &filename[..(filename.len()-4)]).as_str(),])
                    .spawn() {
                        Ok(_) => {},
                        Err(_) => println!("Unable to generate png file (only svg generated). Is 'rsvg-convert' installed?")
                    };
            },
            "" => {},
            _  => { println!("Render format not recognized. Only svg was generated."); }
        }

        Ok(())
    }

    #[inline]
    fn to_y(&self, f: f64) -> f64 {
        let scale = (self.size.1 - self.margins.1.0 - self.margins.1.1) / (self.ylims.1 - self.ylims.0);
        self.size.1 - self.margins.1.0 - (f - self.ylims.0) * scale
    }

    #[inline]
    fn convert_height(&self, h: f64) -> f64 {
        let scale = (self.size.1 - self.margins.1.0 - self.margins.1.1) / (self.ylims.1 - self.ylims.0);
        h * scale
    }

    fn generate_ticks(&self) -> Vec<(f64, String, f64)> {
        match &self.ticks {
            Ticks::LabeledList(tick_list) => {
                tick_list.iter().map(|(val, label)| {
                    (self.to_y(*val), label.clone(), *val)
                }).collect()
            },
            Ticks::Number(nticks) => {
                let delta = (self.ylims.1 - self.ylims.0) / (*nticks - 1) as f64;
                (0..*nticks).map(|i| {
                    let val = self.ylims.0 + (i as f64) * delta;
                    (self.to_y(val), format!("{0:.1$}", val, self.tick_precision), val)
                }).collect()
            },
            Ticks::List(lst) => {
                (*lst).iter().map(|val| {
                    (self.to_y(*val), format!("{0:.1$}", val, self.tick_precision), *val)
                }).collect()
            },
        }
    }

    fn generate_x_positions(&self) -> Vec<f64> {
        let (lw, sw, fw) = (30., 2., 3.);
        let scale = (self.size.0 - self.margins.0.0 - self.margins.0.1) /
                       (  (L as f64 + 2.) * lw
                        + (2_f64.powf(L as f64) - L as f64) * sw
                        +  2_f64.powf(L as f64) * fw             );

        let mut last_g_mutations: i64 = 0;
        let mut x_positions = Vec::with_capacity(self.genotypes.len());
        let mut current_position = self.margins.0.0 + (fw / 2. + lw - sw) * scale;

        for i in 0..self.genotypes.len() {
            if last_g_mutations != self.genotypes[i].sum() as i64 {
                last_g_mutations += 1;
                current_position += lw * scale;
            } else {
                current_position += sw * scale;
            }
            x_positions.push(current_position);
            current_position += fw * scale;
        }
        x_positions
    }

    fn generate_block_positions(&self) -> Vec<(f64, f64)> {
        let (lw, sw, fw) = (30., 2., 3.);
        let scale = (self.size.0 - self.margins.0.0 - self.margins.0.1) /
                       (  (L as f64 + 2.) * lw
                        + (2_f64.powf(L as f64) - L as f64) * sw
                        +  2_f64.powf(L as f64) * fw             );

        let mut last_g_mutations: i64 = 0;
        let mut x_positions = Vec::with_capacity(self.genotypes.len());
        let mut current_position = self.margins.0.0 + (fw / 2. + lw - sw) * scale;

        for i in 0..self.genotypes.len() {
            if last_g_mutations != self.genotypes[i].sum() as i64 {
                last_g_mutations += 1;

                let n = binomial_coefficient((last_g_mutations - 1) as usize, L) as f64;

                x_positions.push((
                    current_position - (lw / 3.0 + n * fw + (n-1.) * sw) * scale,
                    current_position + (lw / 3.0 - fw) * scale
                ));
                current_position += lw * scale;
            } else {
                current_position += sw * scale;
            }
            current_position += fw * scale;
        }

        x_positions.push((
            current_position - (lw / 3.0 + fw) * scale,
            current_position + (lw / 3.0 - fw) * scale
        ));

        x_positions
    }

}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}
fn binomial_coefficient(k: usize, n: usize) -> usize {
    factorial(n) / (factorial(k) * factorial(n - k))
}
